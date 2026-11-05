use nalgebra_sparse::{CscMatrix, CsrMatrix};

use crate::fields::Kind;
use crate::output::LayerMatrix;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Cell {
    Int(Option<i32>),
    Float(Option<f32>),
}

pub struct Matrices {
    pub x: CsrMatrix<i8>,
    pub layers: Vec<(String, LayerMatrix)>,
    pub missing: Vec<(String, u64)>,
}

struct Columns<T> {
    indptr: Vec<usize>,
    indices: Vec<usize>,
    data: Vec<T>,
}

impl<T: nalgebra_sparse::na::Scalar> Columns<T> {
    fn new() -> Self {
        Self {
            indptr: vec![0],
            indices: Vec::new(),
            data: Vec::new(),
        }
    }

    fn push(&mut self, row: u32, value: T) {
        self.indices.push(row as usize);
        self.data.push(value);
    }

    fn end_column(&mut self) {
        self.indptr.push(self.data.len());
    }

    fn into_csr(self, n_obs: usize) -> CsrMatrix<T> {
        let n_vars = self.indptr.len() - 1;
        let csc = CscMatrix::try_from_csc_data(n_obs, n_vars, self.indptr, self.indices, self.data)
            .expect("rows must be sorted, unique and below n_obs within each column");
        CsrMatrix::from(&csc)
    }
}

enum LayerColumns {
    Int(Columns<i32>),
    Float(Columns<f32>),
}

struct Layer {
    name: String,
    columns: LayerColumns,
    missing: u64,
}

pub struct SparseBuilder {
    n_obs: usize,
    x: Columns<i8>,
    layers: Vec<Layer>,
}

impl SparseBuilder {
    /// # Panics
    /// Panics if a layer is neither `Kind::Integer` nor `Kind::Float`.
    #[must_use]
    pub fn new(n_obs: usize, layers: &[(String, Kind)]) -> Self {
        let layers = layers
            .iter()
            .map(|(name, kind)| Layer {
                name: name.clone(),
                columns: match kind {
                    Kind::Integer => LayerColumns::Int(Columns::new()),
                    Kind::Float => LayerColumns::Float(Columns::new()),
                    Kind::Flag | Kind::Text => panic!("layer {name} must be integer or float"),
                },
                missing: 0,
            })
            .collect();
        Self {
            n_obs,
            x: Columns::new(),
            layers,
        }
    }

    /// Rows must ascend within a column.
    pub fn push_copies(&mut self, row: u32, copies: i8) {
        if copies > 0 {
            self.x.push(row, copies);
        }
    }

    /// Rows must ascend within a column. A missing value counts only at a cell with reads.
    ///
    /// # Panics
    /// Panics if a cell's kind differs from its layer's.
    pub fn push_cell(&mut self, layer: usize, row: u32, cell: Cell, has_reads: bool) {
        let layer = &mut self.layers[layer];
        match (&mut layer.columns, cell) {
            (LayerColumns::Int(c), Cell::Int(Some(value))) => {
                if value != 0 {
                    c.push(row, value);
                }
            }
            (LayerColumns::Float(c), Cell::Float(Some(value))) => {
                if value != 0.0 {
                    c.push(row, value);
                }
            }
            (LayerColumns::Int(_), Cell::Int(None))
            | (LayerColumns::Float(_), Cell::Float(None)) => {
                layer.missing += u64::from(has_reads);
            }
            _ => panic!("cell kind does not match layer {}", layer.name),
        }
    }

    pub fn end_column(&mut self) {
        self.x.end_column();
        for layer in &mut self.layers {
            match &mut layer.columns {
                LayerColumns::Int(c) => c.end_column(),
                LayerColumns::Float(c) => c.end_column(),
            }
        }
    }

    #[must_use]
    pub fn finish(self) -> Matrices {
        let n_obs = self.n_obs;
        let mut layers = Vec::with_capacity(self.layers.len());
        let mut missing = Vec::new();
        for layer in self.layers {
            if layer.missing > 0 {
                missing.push((layer.name.clone(), layer.missing));
            }
            let matrix = match layer.columns {
                LayerColumns::Int(c) => LayerMatrix::I32(c.into_csr(n_obs)),
                LayerColumns::Float(c) => LayerMatrix::F32(c.into_csr(n_obs)),
            };
            layers.push((layer.name, matrix));
        }
        Matrices {
            x: self.x.into_csr(n_obs),
            layers,
            missing,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn int_layers() -> Vec<(String, Kind)> {
        vec![("AD_alt".to_string(), Kind::Integer)]
    }

    fn csr_entry<T: Copy>(m: &CsrMatrix<T>, row: usize, col: usize) -> Option<T> {
        let r = m.row(row);
        let at = r.col_indices().iter().position(|&c| c == col)?;
        Some(r.values()[at])
    }

    fn int_cell(value: i32) -> Cell {
        Cell::Int(Some(value))
    }

    #[test]
    fn sparse_builder_stores_only_nonzero_counts() {
        let mut b = SparseBuilder::new(3, &int_layers());
        b.push_cell(0, 0, int_cell(0), true);
        b.push_cell(0, 2, int_cell(4), true);
        b.end_column();
        let m = b.finish();
        let LayerMatrix::I32(ad) = &m.layers[0].1 else {
            panic!("wrong kind")
        };
        assert_eq!(ad.nnz(), 1);
        assert_eq!(csr_entry(ad, 0, 0), None);
        assert_eq!(csr_entry(ad, 2, 0), Some(4));
        assert_eq!(m.missing.len(), 0);
    }

    #[test]
    fn a_zero_float_is_not_stored() {
        let mut b = SparseBuilder::new(2, &[("AF".to_string(), Kind::Float)]);
        b.push_cell(0, 0, Cell::Float(Some(0.0)), true);
        b.push_cell(0, 1, Cell::Float(Some(0.5)), true);
        b.end_column();
        let m = b.finish();
        let LayerMatrix::F32(af) = &m.layers[0].1 else {
            panic!("wrong kind")
        };
        assert_eq!(af.nnz(), 1);
        assert_eq!(csr_entry(af, 1, 0), Some(0.5));
    }

    #[test]
    fn x_stores_only_copies_above_zero() {
        let mut b = SparseBuilder::new(3, &int_layers());
        b.push_copies(0, 0);
        b.push_copies(1, 1);
        b.push_copies(2, 2);
        b.end_column();
        let m = b.finish();
        assert_eq!(m.x.nnz(), 2);
        assert_eq!(csr_entry(&m.x, 0, 0), None);
        assert_eq!(csr_entry(&m.x, 1, 0), Some(1));
        assert_eq!(csr_entry(&m.x, 2, 0), Some(2));
    }

    #[test]
    fn x_does_not_follow_the_layers() {
        let mut b = SparseBuilder::new(2, &int_layers());
        b.push_copies(0, 1);
        b.push_cell(0, 1, int_cell(3), true);
        b.end_column();
        let m = b.finish();
        assert_eq!(csr_entry(&m.x, 0, 0), Some(1));
        let LayerMatrix::I32(ad) = &m.layers[0].1 else {
            panic!("wrong kind")
        };
        assert_eq!(csr_entry(ad, 0, 0), None);
        assert_eq!(csr_entry(ad, 1, 0), Some(3));
    }

    #[test]
    #[should_panic(expected = "cell kind does not match layer AD_alt")]
    fn a_cell_must_match_its_layer() {
        let mut b = SparseBuilder::new(1, &int_layers());
        b.push_cell(0, 0, Cell::Float(Some(1.0)), true);
    }

    #[test]
    fn a_missing_int_is_unstored_and_counted_at_a_cell_with_reads() {
        let mut b = SparseBuilder::new(2, &int_layers());
        b.push_cell(0, 0, Cell::Int(None), true);
        b.push_cell(0, 1, int_cell(3), true);
        b.end_column();
        let m = b.finish();
        let LayerMatrix::I32(ad) = &m.layers[0].1 else {
            panic!("wrong kind")
        };
        assert_eq!(ad.nnz(), 1);
        assert_eq!(m.missing, vec![("AD_alt".to_string(), 1)]);
    }

    #[test]
    fn a_missing_float_is_unstored_and_counted_at_a_cell_with_reads() {
        let mut b = SparseBuilder::new(1, &[("AF".to_string(), Kind::Float)]);
        b.push_cell(0, 0, Cell::Float(None), true);
        b.end_column();
        let m = b.finish();
        let LayerMatrix::F32(af) = &m.layers[0].1 else {
            panic!("wrong kind")
        };
        assert_eq!(af.nnz(), 0);
        assert_eq!(m.missing, vec![("AF".to_string(), 1)]);
    }

    #[test]
    fn a_missing_value_at_a_cell_with_no_reads_is_not_counted() {
        let mut b = SparseBuilder::new(1, &int_layers());
        b.push_cell(0, 0, Cell::Int(None), false);
        b.end_column();
        assert_eq!(b.finish().missing.len(), 0);
    }

    #[test]
    fn columns_accumulate_across_pushes() {
        let mut b = SparseBuilder::new(2, &int_layers());
        b.push_copies(1, 1);
        b.push_cell(0, 1, int_cell(7), true);
        b.end_column();
        b.push_copies(0, 3);
        b.push_copies(1, 2);
        b.push_cell(0, 0, int_cell(8), true);
        b.push_cell(0, 1, int_cell(9), true);
        b.end_column();
        let m = b.finish();
        assert_eq!((m.x.nrows(), m.x.ncols()), (2, 2));
        assert_eq!(csr_entry(&m.x, 1, 0), Some(1));
        assert_eq!(csr_entry(&m.x, 0, 1), Some(3));
        let LayerMatrix::I32(ad) = &m.layers[0].1 else {
            panic!("wrong kind")
        };
        assert_eq!(csr_entry(ad, 1, 1), Some(9));
    }
}
