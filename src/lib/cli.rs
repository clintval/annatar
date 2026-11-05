//! Argument parsing, help and usage errors for the `annatar` command.
use std::ffi::{OsStr, OsString};
use std::io::{IsTerminal, Write};
use std::num::NonZero;
use std::path::PathBuf;

use clap::builder::styling::{AnsiColor, Effects, Style, Styles};
use clap::builder::{OsStringValueParser, TypedValueParser};
use clap::{Args, CommandFactory, FromArgMatches, Parser, Subcommand};
use env_logger::{Env, Logger};

use crate::convert::Config;
use crate::fields::{Exclusion, Exclusions};

const HEADER: Style = AnsiColor::Green.on_default().effects(Effects::BOLD);
const USAGE: Style = AnsiColor::Green.on_default().effects(Effects::BOLD);
const LITERAL: Style = AnsiColor::Cyan.on_default().effects(Effects::BOLD);
const PLACEHOLDER: Style = AnsiColor::Cyan.on_default();
const ERROR: Style = AnsiColor::Red.on_default().effects(Effects::BOLD);
const VALID: Style = AnsiColor::Cyan.on_default().effects(Effects::BOLD);
const INVALID: Style = AnsiColor::Yellow.on_default().effects(Effects::BOLD);
const TITLE: Style = AnsiColor::Magenta.on_default();
const CODE: Style = AnsiColor::Yellow.on_default();
const FOOTER: Style = AnsiColor::Magenta.on_default();
const FADED: Style = AnsiColor::BrightBlack.on_default();

const CARGO_STYLING: Styles = Styles::styled()
    .header(HEADER)
    .usage(USAGE)
    .literal(LITERAL)
    .placeholder(PLACEHOLDER)
    .error(ERROR)
    .valid(VALID)
    .invalid(INVALID);

/// Convert a multi-sample VCF of single cells into an h5ad.
///
/// Each VCF sample becomes a cell in `obs` and each ALT allele a variant in
/// `var`. The `X` matrix holds the
/// number of copies of this ALT in the called genotype; hom-ref and uncalled
/// are not stored. Numeric FORMAT fields become sparse
/// layers of plain counts, with zeros and missing values not stored.
/// VEP annotation is optional: with `INFO/CSQ`, one consequence is picked per
/// ALT and becomes typed `var` columns. The file is written once the whole VCF
/// has been read, so a failed run leaves no output.
///
/// EXAMPLES
///
///   annatar 'cells.vcf.gz' 'cells.h5ad'
///   annatar --genes 'panel-genes.txt' 'cells.vep.vcf.gz' 'cells.h5ad'
///   bcftools view 'cells.bcf' | annatar - 'cells.h5ad'
///   annatar merge --output 'merged.h5ad' 'a.h5ad' 'b.h5ad'
#[derive(Debug, Parser)]
#[command(
    author,
    version,
    rename_all = "kebab-case",
    verbatim_doc_comment,
    arg_required_else_help = true,
    args_conflicts_with_subcommands = true,
    subcommand_negates_reqs = true,
    disable_help_subcommand = true
)]
#[clap(styles = CARGO_STYLING)]
struct Cli {
    #[command(flatten)]
    convert: ConvertArgs,

    #[command(subcommand)]
    command: Option<Command>,
}

#[derive(Debug, Subcommand)]
enum Command {
    Merge(MergeArgs),
}

/// Merge annatar's h5ad files, one per library, into one h5ad.
///
/// Each input is one library, named after its file. Barcodes gain the library
/// name, `obs["library"]` names the library, and variants are the union, sorted
/// by VCF header contig order, `pos`, `ref` and `alt`. A shared variant takes its
/// whole `var` row from the first library holding it, and `n_cells_covered`,
/// `n_cells_alt`, `ref_reads` and `alt_reads` are recounted from the merged
/// layers. The file is written once every input has been read, so a failed run
/// leaves no output.
///
/// EXAMPLES
///
///   annatar merge --output 'merged.h5ad' 'a.h5ad' 'b.h5ad'
#[derive(Debug, Args)]
#[command(
    rename_all = "kebab-case",
    verbatim_doc_comment,
    arg_required_else_help = true
)]
struct MergeArgs {
    /// Output h5ad. HDF5 cannot stream, so this must be a path.
    #[arg(
        short,
        long,
        value_name = "H5AD",
        value_parser = OsStringValueParser::new().try_map(output_path),
        verbatim_doc_comment
    )]
    output: PathBuf,

    /// Input h5ad files written by annatar, one per library.
    #[arg(value_name = "H5AD", required = true, verbatim_doc_comment)]
    inputs: Vec<PathBuf>,
}

#[derive(Debug, Args)]
struct ConvertArgs {
    /// Input VCF or BCF, plain or bgzipped, or `-` for standard input.
    ///
    /// A file named `merge` must be written `./merge`.
    #[arg(value_name = "VCF", required = true, verbatim_doc_comment)]
    input: Option<PathBuf>,

    /// Output h5ad. HDF5 cannot stream, so this must be a path.
    #[arg(
        value_name = "H5AD",
        required = true,
        value_parser = OsStringValueParser::new().try_map(output_path),
        verbatim_doc_comment
    )]
    output: Option<PathBuf>,

    /// Panel genes preferred by the pick, one HGNC symbol or Ensembl gene ID per
    /// line.
    #[arg(long, value_name = "FILE", verbatim_doc_comment)]
    genes: Option<PathBuf>,

    /// NCBI MANE summary, plain or gzipped, used when CSQ has no MANE subfields.
    #[arg(long, value_name = "FILE", verbatim_doc_comment)]
    mane_table: Option<PathBuf>,

    /// Field to leave out, named `INFO/ID`, `FORMAT/ID` or `INFO/CSQ/subfield`.
    ///
    /// Repeat the flag to leave out several, such as
    /// `INFO/CSQ/Existing_variation` and `INFO/CSQ/SOMATIC`. A name must match
    /// the header exactly: one it does not declare is an error, except that
    /// `INFO/CSQ/subfield` only warns when there is no `INFO/CSQ`. Excluding
    /// `INFO/CSQ` drops all annotation.
    #[arg(long, value_name = "FIELD", verbatim_doc_comment)]
    exclude_field: Vec<Exclusion>,

    /// BGZF decompression threads.
    #[arg(
        short = '@',
        long,
        value_name = "N",
        default_value = "1",
        verbatim_doc_comment
    )]
    threads: NonZero<usize>,
}

fn output_path(path: OsString) -> Result<PathBuf, String> {
    if path == "-" {
        Err("HDF5 cannot stream, so the output must be a file, not -".to_string())
    } else {
        Ok(PathBuf::from(path))
    }
}

fn wants_color(no_color: Option<&OsStr>, is_terminal: bool) -> bool {
    is_terminal && no_color.is_none_or(OsStr::is_empty)
}

fn rendered(error: &clap::Error, color: bool) -> String {
    let styled = error.render();
    if color {
        styled.ansi().to_string()
    } else {
        styled.to_string()
    }
}

// Writing directly ignores a closed pipe, which print! would turn into a panic.
fn report(error: &clap::Error, stdout_color: bool, stderr_color: bool) -> i32 {
    let (mut stream, color): (Box<dyn Write>, bool) = if error.use_stderr() {
        (Box::new(std::io::stderr()), stderr_color)
    } else {
        (Box::new(std::io::stdout()), stdout_color)
    };
    let _ = stream.write_all(rendered(error, color).as_bytes());
    let _ = stream.flush();
    error.exit_code()
}

fn title_first_line(text: &str) -> String {
    let title = TITLE.render();
    let reset = TITLE.render_reset();
    match text.split_once('\n') {
        Some((first, rest)) => format!("{title}{first}{reset}\n{rest}"),
        None => format!("{title}{text}{reset}"),
    }
}

fn paint_backtick_terms(text: &str, code: &str, after: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut rest = text;
    while let Some(open) = rest.find('`') {
        out.push_str(&rest[..open]);
        let tail = &rest[open + 1..];
        if let Some(close) = tail.find('`') {
            out.push_str(code);
            out.push_str(&tail[..close]);
            out.push_str(after);
            rest = &tail[close + 1..];
        } else {
            out.push('`');
            rest = tail;
            break;
        }
    }
    out.push_str(rest);
    out
}

fn split_two_column(content: &str) -> Option<(&str, &str, &str)> {
    let bytes = content.as_bytes();
    let mut seen_term = false;
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b' ' {
            let start = i;
            while i < bytes.len() && bytes[i] == b' ' {
                i += 1;
            }
            if seen_term && i - start >= 2 {
                return Some((&content[..start], &content[start..i], &content[i..]));
            }
        } else {
            seen_term = true;
            i += 1;
        }
    }
    None
}

fn style_help_text(text: &str) -> String {
    let code = CODE.render().to_string();
    let faded = FADED.render().to_string();
    let reset = CODE.render_reset().to_string();
    text.split_inclusive('\n')
        .map(|line| {
            let (content, newline) = match line.strip_suffix('\n') {
                Some(content) => (content, "\n"),
                None => (line, ""),
            };
            if !content.starts_with("  ") || content.trim().is_empty() {
                return format!("{}{newline}", paint_backtick_terms(content, &code, &reset));
            }
            let indent = content.len() - content.trim_start().len();
            if indent > 2 {
                let painted = paint_backtick_terms(content, &code, &faded);
                format!("{faded}{painted}{reset}{newline}")
            } else if content.trim_start().starts_with(env!("CARGO_PKG_NAME")) {
                let painted = paint_backtick_terms(content, &code, &code);
                format!("{code}{painted}{reset}{newline}")
            } else if let Some((left, gap, right)) = split_two_column(content) {
                let left = paint_backtick_terms(left, &code, &code);
                let right = paint_backtick_terms(right, &code, &faded);
                format!("{code}{left}{reset}{gap}{faded}{right}{reset}{newline}")
            } else {
                let painted = paint_backtick_terms(content, &code, &code);
                format!("{code}{painted}{reset}{newline}")
            }
        })
        .collect()
}

fn decorate_help(cmd: clap::Command) -> clap::Command {
    let about = cmd.get_about().map(ToString::to_string);
    let long_about = cmd
        .get_long_about()
        .map(ToString::to_string)
        .or_else(|| about.clone());

    let footer = format!(
        "{}MIT License 2026 · Clint Valentine{}",
        FOOTER.render(),
        FOOTER.render_reset(),
    );

    let subcommands: Vec<String> = cmd
        .get_subcommands()
        .map(|subcommand| subcommand.get_name().to_string())
        .collect();
    let mut cmd = subcommands
        .into_iter()
        .fold(cmd, |cmd, name| cmd.mut_subcommand(name, decorate_help));
    if let Some(about) = about {
        cmd = cmd.about(title_first_line(&style_help_text(&about)));
    }
    if let Some(long_about) = long_about {
        cmd = cmd.long_about(title_first_line(&style_help_text(&long_about)));
    }
    cmd = cmd
        // Help text is hand-wrapped, and clap would count the injected ANSI escapes as width.
        .term_width(usize::MAX)
        .next_line_help(true)
        .after_help(footer.clone())
        .after_long_help(footer);
    cmd.mut_args(|mut arg| {
        if let Some(help) = arg.get_help().map(ToString::to_string) {
            arg = arg.help(style_help_text(&help));
        }
        if let Some(long_help) = arg.get_long_help().map(ToString::to_string) {
            arg = arg.long_help(style_help_text(&long_help));
        }
        arg
    })
}

pub enum Invocation {
    Convert {
        config: Config,
        output: PathBuf,
    },
    Merge {
        inputs: Vec<PathBuf>,
        output: PathBuf,
    },
}

pub struct Parsed {
    pub invocation: Invocation,
    pub logger: Logger,
}

fn stderr_logger(color: bool) -> Logger {
    let write_style = if color {
        env_logger::WriteStyle::Always
    } else {
        env_logger::WriteStyle::Never
    };
    env_logger::Builder::from_env(Env::default().default_filter_or("info"))
        .write_style(write_style)
        .build()
}

/// Parses `args`, or renders help or a usage error and returns the exit code.
pub fn parse(args: &[OsString]) -> Result<Parsed, i32> {
    let no_color = std::env::var_os("NO_COLOR");
    let stdout_color = wants_color(no_color.as_deref(), std::io::stdout().is_terminal());
    let stderr_color = wants_color(no_color.as_deref(), std::io::stderr().is_terminal());

    let mut cmd = decorate_help(Cli::command());
    let cli = cmd
        .try_get_matches_from_mut(args)
        .and_then(|matches| Cli::from_arg_matches(&matches))
        .map_err(|error| report(&error, stdout_color, stderr_color))?;

    let invocation = if let Some(Command::Merge(args)) = cli.command {
        Invocation::Merge {
            inputs: args.inputs,
            output: args.output,
        }
    } else {
        let convert = cli.convert;
        Invocation::Convert {
            config: Config {
                input: convert.input.expect("clap requires the input"),
                genes: convert.genes,
                mane_table: convert.mane_table,
                exclusions: Exclusions::new(convert.exclude_field),
                threads: convert.threads,
                command_line: args
                    .iter()
                    .map(|arg| arg.to_string_lossy().into_owned())
                    .collect::<Vec<_>>()
                    .join(" "),
            },
            output: convert.output.expect("clap requires the output"),
        }
    };
    Ok(Parsed {
        invocation,
        logger: stderr_logger(stderr_color),
    })
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[test]
    fn the_output_must_be_a_file() {
        assert!(output_path("-".into()).is_err());
        assert_eq!(output_path("out.h5ad".into()), Ok("out.h5ad".into()));
    }

    #[test]
    fn two_column_rows_split_at_the_first_wide_gap() {
        assert_eq!(
            split_two_column("  1           no copies"),
            Some(("  1", "           ", "no copies"))
        );
        assert_eq!(split_two_column("  annatar -o out.h5ad"), None);
    }

    #[rstest]
    #[case(None, true, true)]
    #[case(None, false, false)]
    #[case(Some(""), true, true)]
    #[case(Some(""), false, false)]
    #[case(Some("1"), true, false)]
    #[case(Some("1"), false, false)]
    #[case(Some("0"), true, false)]
    #[case(Some("0"), false, false)]
    fn color_needs_a_terminal_and_no_non_empty_no_color(
        #[case] no_color: Option<&str>,
        #[case] is_terminal: bool,
        #[case] expected: bool,
    ) {
        assert_eq!(wants_color(no_color.map(OsStr::new), is_terminal), expected);
    }

    #[test]
    fn a_rendered_error_has_escapes_only_when_asked() {
        let error = decorate_help(Cli::command())
            .try_get_matches_from(["annatar", "--bogus"])
            .unwrap_err();
        assert!(rendered(&error, true).contains('\u{1b}'));
        assert!(!rendered(&error, false).contains('\u{1b}'));
    }

    #[test]
    fn rendered_help_loses_its_painted_escapes_without_color() {
        let error = decorate_help(Cli::command())
            .try_get_matches_from(["annatar", "--help"])
            .unwrap_err();
        assert!(rendered(&error, true).contains(&*TITLE.render().to_string()));
        let plain = rendered(&error, false);
        assert!(!plain.contains('\u{1b}'));
        assert!(plain.contains("MIT License"));
    }

    #[test]
    fn help_text_drops_backticks_and_paints_terms() {
        let text = "Use `INFO/CSQ/SYMBOL`.\n  1  one copy\n";
        let painted = style_help_text(text);
        assert!(!painted.contains('`'));
        assert!(painted.contains(&format!("{}INFO/CSQ/SYMBOL", CODE.render())));
    }
}
