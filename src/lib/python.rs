//! The `annatar._native` extension module: a conversion handed to Python as `NumPy` arrays.

use std::cell::RefCell;
use std::collections::HashMap;
use std::ffi::{CStr, CString, OsString};
use std::io;
use std::num::NonZeroUsize;
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

use env_logger::Logger;
use log::{Level, LevelFilter, Log, Metadata, Record, error, info};
use nalgebra_sparse::CsrMatrix;
use nalgebra_sparse::na::Scalar;
use numpy::{Element, IntoPyArray, PyArray1};
use pyo3::IntoPyObjectExt;
use pyo3::exceptions::{
    PyFileNotFoundError, PyOSError, PyRuntimeWarning, PyUserWarning, PyValueError,
};
use pyo3::prelude::*;

use crate::cli::{Invocation, Parsed, parse};
use crate::convert::{Config, convert};
use crate::fields::{Exclusion, Exclusions};
use crate::order::variant_order;
use crate::output::{Column, LayerMatrix, Output, UnsValue};

type Sparse<'py> = (Bound<'py, PyAny>, Bound<'py, PyAny>, Bound<'py, PyAny>);
type NamedColumn<'py> = (String, &'static str, Bound<'py, PyAny>, Bound<'py, PyAny>);

thread_local! {
    static MESSAGES: RefCell<Vec<(Level, String)>> = const { RefCell::new(Vec::new()) };
}

const LOGGER_TAKEN: &CStr = c"another Rust logger is installed, so annatar warnings are lost";

struct MessageLog;

static MESSAGE_LOG: MessageLog = MessageLog;
static STDERR_LOG: OnceLock<Logger> = OnceLock::new();

impl Log for MessageLog {
    fn enabled(&self, metadata: &Metadata) -> bool {
        match STDERR_LOG.get() {
            Some(logger) => logger.enabled(metadata),
            None => metadata.level() <= Level::Info && metadata.target().starts_with("annatar"),
        }
    }

    fn log(&self, record: &Record) {
        if let Some(logger) = STDERR_LOG.get() {
            logger.log(record);
        } else if self.enabled(record.metadata()) {
            MESSAGES.with_borrow_mut(|messages| {
                messages.push((record.level(), record.args().to_string()));
            });
        }
    }

    fn flush(&self) {}
}

// A second import of this module runs its setup again, and that must keep the first answer.
fn install_logger() -> bool {
    static INSTALLED: OnceLock<bool> = OnceLock::new();
    *INSTALLED.get_or_init(|| {
        let installed = log::set_logger(&MESSAGE_LOG).is_ok();
        if installed {
            log::set_max_level(LevelFilter::Info);
        }
        installed
    })
}

fn to_python(error: &anyhow::Error) -> PyErr {
    let message = format!("{error:#}");
    let kind = error
        .chain()
        .find_map(|cause| cause.downcast_ref::<io::Error>())
        .map(io::Error::kind);
    match kind {
        Some(io::ErrorKind::NotFound) => PyFileNotFoundError::new_err(message),
        None | Some(io::ErrorKind::InvalidData | io::ErrorKind::UnexpectedEof) => {
            PyValueError::new_err(message)
        }
        Some(_) => PyOSError::new_err(message),
    }
}

// Stack level 2 is the caller of the Python function that wraps this module.
fn emit_messages(py: Python<'_>) -> PyResult<()> {
    let category = py.get_type::<PyUserWarning>();
    let logger = py
        .import("logging")?
        .call_method1("getLogger", ("annatar",))?;
    for (level, message) in MESSAGES.take() {
        if level <= Level::Warn {
            PyErr::warn(py, &category, &CString::new(message)?, 2)?;
        } else {
            logger.call_method1("info", (message,))?;
        }
    }
    Ok(())
}

// scipy copies index arrays wider than a matrix needs, so these arrive at the narrowest width.
fn index_array(py: Python<'_>, values: Vec<usize>) -> PyResult<Bound<'_, PyAny>> {
    let narrow: Result<Vec<i32>, _> = values.iter().map(|&v| i32::try_from(v)).collect();
    if let Ok(narrow) = narrow {
        narrow.into_pyarray(py).into_bound_py_any(py)
    } else {
        let wide: Vec<isize> = values.into_iter().map(usize::cast_signed).collect();
        wide.into_pyarray(py).into_bound_py_any(py)
    }
}

fn sparse<T: Element + Scalar>(py: Python<'_>, matrix: CsrMatrix<T>) -> PyResult<Sparse<'_>> {
    let (indptr, indices, data) = matrix.disassemble();
    Ok((
        data.into_pyarray(py).into_any(),
        index_array(py, indices)?,
        index_array(py, indptr)?,
    ))
}

fn nullable<'py, T: Copy + Default + Element>(
    py: Python<'py>,
    values: &[Option<T>],
) -> (Bound<'py, PyAny>, Bound<'py, PyAny>) {
    let filled: Vec<T> = values.iter().map(|v| v.unwrap_or_default()).collect();
    let mask: Vec<bool> = values.iter().map(Option::is_none).collect();
    (
        filled.into_pyarray(py).into_any(),
        mask.into_pyarray(py).into_any(),
    )
}

// Categories keep the order of first appearance, as the h5ad writer's do.
fn categorical(values: Vec<Option<String>>) -> PyResult<(Vec<i32>, Vec<String>)> {
    let mut codes = HashMap::new();
    let mut categories = Vec::new();
    let codes = values
        .into_iter()
        .map(|value| {
            let Some(value) = value else { return Ok(-1) };
            if let Some(code) = codes.get(&value) {
                return Ok(*code);
            }
            let code = i32::try_from(categories.len())?;
            codes.insert(value.clone(), code);
            categories.push(value);
            Ok(code)
        })
        .collect::<Result<_, std::num::TryFromIntError>>()
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
    Ok((codes, categories))
}

fn named_column(py: Python<'_>, (name, column): (String, Column)) -> PyResult<NamedColumn<'_>> {
    let none = || py.None().into_bound(py);
    Ok(match column {
        Column::Int(values) => (name, "int64", values.into_pyarray(py).into_any(), none()),
        Column::Float(values) => (name, "float64", values.into_pyarray(py).into_any(), none()),
        Column::Flag(values) => (name, "bool", values.into_pyarray(py).into_any(), none()),
        Column::NullableInt(values) => {
            let (values, mask) = nullable(py, &values);
            (name, "Int64", values, mask)
        }
        Column::NullableFlag(values) => {
            let (values, mask) = nullable(py, &values);
            (name, "boolean", values, mask)
        }
        Column::Text(values) => {
            let (codes, categories) = categorical(values)?;
            let categories = categories.into_bound_py_any(py)?;
            (
                name,
                "category",
                codes.into_pyarray(py).into_any(),
                categories,
            )
        }
    })
}

fn columns(py: Python<'_>, columns: Vec<(String, Column)>) -> PyResult<Vec<NamedColumn<'_>>> {
    columns
        .into_iter()
        .map(|column| named_column(py, column))
        .collect()
}

#[derive(IntoPyObject)]
struct Conversion<'py> {
    obs_names: Vec<String>,
    obs: Vec<NamedColumn<'py>>,
    var_names: Vec<String>,
    var: Vec<NamedColumn<'py>>,
    x: Sparse<'py>,
    layers: Vec<(String, Sparse<'py>)>,
    uns: Vec<(String, Bound<'py, PyAny>)>,
}

impl<'py> Conversion<'py> {
    fn new(py: Python<'py>, output: Output) -> PyResult<Self> {
        let layers = output
            .layers
            .into_iter()
            .map(|(name, matrix)| {
                let arrays = match matrix {
                    LayerMatrix::I32(matrix) => sparse(py, matrix)?,
                    LayerMatrix::F32(matrix) => sparse(py, matrix)?,
                };
                Ok((name, arrays))
            })
            .collect::<PyResult<_>>()?;
        let uns = output
            .uns
            .into_iter()
            .map(|(name, value)| {
                let value = match value {
                    UnsValue::Text(text) => text.into_bound_py_any(py)?,
                    UnsValue::List(values) => values.into_bound_py_any(py)?,
                };
                Ok((name, value))
            })
            .collect::<PyResult<_>>()?;
        Ok(Self {
            obs_names: output.obs_names,
            obs: columns(py, output.obs)?,
            var_names: output.var_names,
            var: columns(py, output.var)?,
            x: sparse(py, output.x)?,
            layers,
            uns,
        })
    }
}

// PyO3 extracts arguments by value.
#[allow(clippy::needless_pass_by_value)]
#[pyfunction(name = "variant_order")]
fn order_variants(
    py: Python<'_>,
    ranks: Vec<i64>,
    positions: Vec<i64>,
    references: Vec<String>,
    alternates: Vec<String>,
) -> PyResult<Bound<'_, PyArray1<usize>>> {
    let order = py.detach(|| variant_order(&ranks, &positions, &references, &alternates));
    let order = order.map_err(|error| PyOSError::new_err(error.to_string()))?;
    Ok(order.into_pyarray(py))
}

#[pyfunction]
#[pyo3(signature = (path, *, genes, mane_table, exclude_fields, threads, command_line))]
fn read_vcf(
    py: Python<'_>,
    path: PathBuf,
    genes: Option<PathBuf>,
    mane_table: Option<PathBuf>,
    exclude_fields: Vec<String>,
    threads: NonZeroUsize,
    command_line: String,
) -> PyResult<Conversion<'_>> {
    let exclusions = exclude_fields
        .into_iter()
        .map(|field| field.parse())
        .collect::<Result<Vec<Exclusion>, String>>()
        .map_err(PyValueError::new_err)?;
    let config = Config {
        input: path,
        genes,
        mane_table,
        exclusions: Exclusions::new(exclusions),
        threads,
        command_line,
    };
    MESSAGES.take();
    let output = py.detach(|| convert(&config));
    emit_messages(py)?;
    Conversion::new(py, output.map_err(|error| to_python(&error))?)
}

fn finish(py: Python<'_>, outcome: PyResult<Bound<'_, PyAny>>, output: &Path) -> PyResult<i32> {
    match outcome {
        Ok(counts) => {
            let (cells, variants): (usize, usize) = counts.extract()?;
            let output = output.display();
            info!(target: "annatar", "wrote {cells} cells and {variants} variants to {output}");
            Ok(0)
        }
        Err(failure)
            if failure.is_instance_of::<PyValueError>(py)
                || failure.is_instance_of::<PyOSError>(py) =>
        {
            error!(target: "annatar", "{}", failure.value(py));
            Ok(1)
        }
        Err(failure) => Err(failure),
    }
}

// Python only sees Ctrl-C between bytecodes, so the signal's default action covers the conversion.
fn interruptible(py: Python<'_>, config: &Config) -> PyResult<anyhow::Result<Output>> {
    let signal = py.import("signal")?;
    let interrupt = signal.getattr("SIGINT")?;
    let handler = signal.call_method1("signal", (&interrupt, signal.getattr("SIG_DFL")?))?;
    let result = py.detach(|| convert(config));
    signal.call_method1("signal", (interrupt, handler))?;
    Ok(result)
}

fn dispatch(
    py: Python<'_>,
    args: &[OsString],
    write_conversion: &Bound<'_, PyAny>,
    merge: &Bound<'_, PyAny>,
) -> PyResult<i32> {
    let Parsed { invocation, logger } = match parse(args) {
        Ok(parsed) => parsed,
        Err(code) => return Ok(code),
    };
    log::set_max_level(logger.filter());
    let _ = STDERR_LOG.set(logger);
    match invocation {
        Invocation::Merge { inputs, output } => finish(py, merge.call1((inputs, &output)), &output),
        Invocation::Convert { config, output } => match interruptible(py, &config)? {
            Ok(result) => {
                let conversion = Conversion::new(py, result)?;
                finish(py, write_conversion.call1((conversion, &output)), &output)
            }
            Err(failure) => {
                error!(target: "annatar", "{failure:#}");
                Ok(1)
            }
        },
    }
}

// Rust's panic hook has already printed the message, and 101 is the exit code a Rust binary uses.
#[allow(clippy::needless_pass_by_value)]
#[pyfunction(name = "main")]
fn run_command(
    py: Python<'_>,
    args: Vec<OsString>,
    write_conversion: &Bound<'_, PyAny>,
    merge: &Bound<'_, PyAny>,
) -> PyResult<i32> {
    catch_unwind(AssertUnwindSafe(|| {
        dispatch(py, &args, write_conversion, merge)
    }))
    .unwrap_or(Ok(101))
}

#[pymodule]
fn _native(module: &Bound<'_, PyModule>) -> PyResult<()> {
    if !install_logger() {
        let py = module.py();
        let category = py.get_type::<PyRuntimeWarning>();
        PyErr::warn(py, &category, LOGGER_TAKEN, 1)?;
    }
    module.add_function(wrap_pyfunction!(read_vcf, module)?)?;
    module.add_function(wrap_pyfunction!(run_command, module)?)?;
    module.add_function(wrap_pyfunction!(order_variants, module)?)?;
    module.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
