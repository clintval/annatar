//! Convert a VEP-annotated multi-sample VCF of single cells into an `AnnData` h5ad.
use std::ffi::OsString;
use std::num::NonZero;
use std::path::PathBuf;
use std::process;

use annatar::convert::{Config, convert};
use annatar::fields::{Exclusion, Exclusions};
use anyhow::Error;
use clap::builder::styling::{AnsiColor, Effects, Style, Styles};
use clap::builder::{OsStringValueParser, TypedValueParser};
use clap::{CommandFactory, FromArgMatches, Parser};
use env_logger::Env;
use log::{error, info};

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

/// Convert a VEP-annotated multi-sample VCF of single cells into an h5ad.
///
/// Each VCF sample becomes a cell in `obs` and each ALT allele a variant in
/// `var`. `X` holds a genotype code per cell and variant, numeric FORMAT
/// fields become sparse layers, and one VEP consequence is picked per ALT.
/// The file is written once the whole VCF has been read, so a failed run
/// leaves no output.
///
/// X
///
///   not stored  no genotype
///   1           no copies of this ALT
///   2           one copy of this ALT
///   3           two or more copies of this ALT
///
/// EXAMPLES
///
///   annatar --output 'cells.h5ad' 'cells.vep.vcf.gz'
///   annatar --genes 'panel-genes.txt' --output 'cells.h5ad' 'cells.vep.vcf.gz'
///   bcftools view 'cells.vep.bcf' | annatar --output 'cells.h5ad'
#[derive(Debug, Parser)]
#[command(
    author,
    version,
    color = clap::ColorChoice::Always,
    rename_all = "kebab-case",
    verbatim_doc_comment,
    arg_required_else_help = true
)]
#[clap(styles = CARGO_STYLING)]
struct Cli {
    /// Input VCF or BCF, plain or bgzipped, or `-` for standard input.
    #[arg(value_name = "VCF", default_value = "-", verbatim_doc_comment)]
    input: PathBuf,

    /// Output h5ad. HDF5 cannot stream, so this must be a path.
    #[arg(
        short = 'o',
        long,
        value_name = "H5AD",
        value_parser = OsStringValueParser::new().try_map(output_path),
        verbatim_doc_comment
    )]
    output: PathBuf,

    /// Panel genes preferred by the pick, one HGNC symbol or Ensembl gene ID per
    /// line.
    #[arg(long, value_name = "FILE", verbatim_doc_comment)]
    genes: Option<PathBuf>,

    /// NCBI MANE summary, plain or gzipped, used when CSQ has no MANE subfields.
    #[arg(long, value_name = "FILE", verbatim_doc_comment)]
    mane_table: Option<PathBuf>,

    /// Field to leave out, named `INFO/ID`, `FORMAT/ID` or `CSQ/subfield`.
    ///
    /// Repeat the flag to leave out several, such as
    /// `CSQ/Existing_variation`, `CSQ/SOMATIC` and `CSQ/PHENO`. A name must match
    /// the header exactly: one it does not declare is an error.
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

fn esc(style: Style, color: bool) -> String {
    if color {
        style.render().to_string()
    } else {
        String::new()
    }
}

fn esc_reset(style: Style, color: bool) -> String {
    if color {
        style.render_reset().to_string()
    } else {
        String::new()
    }
}

fn title_first_line(text: &str, color: bool) -> String {
    let title = esc(TITLE, color);
    let reset = esc_reset(TITLE, color);
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

fn style_help_text(text: &str, color: bool) -> String {
    let code = esc(CODE, color);
    let faded = esc(FADED, color);
    let reset = esc_reset(CODE, color);
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

fn decorate_help(cmd: clap::Command, color: bool) -> clap::Command {
    let about = cmd.get_about().map(ToString::to_string);
    let long_about = cmd
        .get_long_about()
        .map(ToString::to_string)
        .or_else(|| about.clone());

    let footer = format!(
        "{}MIT License 2026 · Clint Valentine{}",
        esc(FOOTER, color),
        esc_reset(FOOTER, color),
    );

    let choice = if color {
        clap::ColorChoice::Always
    } else {
        clap::ColorChoice::Never
    };
    let mut cmd = cmd.color(choice);
    if let Some(about) = about {
        cmd = cmd.about(title_first_line(&style_help_text(&about, color), color));
    }
    if let Some(long_about) = long_about {
        cmd = cmd.long_about(title_first_line(
            &style_help_text(&long_about, color),
            color,
        ));
    }
    cmd = cmd
        .next_line_help(true)
        .after_help(footer.clone())
        .after_long_help(footer);
    cmd.mut_args(|mut arg| {
        if let Some(help) = arg.get_help().map(ToString::to_string) {
            arg = arg.help(style_help_text(&help, color));
        }
        if let Some(long_help) = arg.get_long_help().map(ToString::to_string) {
            arg = arg.long_help(style_help_text(&long_help, color));
        }
        arg
    })
}

#[cfg(not(tarpaulin_include))]
fn main() -> Result<(), Error> {
    let color = std::env::var_os("NO_COLOR").is_none();

    let env = Env::default().default_filter_or("info");
    let write_style = if color {
        env_logger::WriteStyle::Auto
    } else {
        env_logger::WriteStyle::Never
    };
    env_logger::Builder::from_env(env)
        .write_style(write_style)
        .init();

    // Help text is hand-wrapped, and clap would count the injected ANSI escapes as width.
    let mut cmd = decorate_help(Cli::command().term_width(usize::MAX), color);
    let matches = cmd.get_matches_mut();
    let cli = Cli::from_arg_matches(&matches).unwrap_or_else(|e| e.exit());

    let config = Config {
        input: cli.input,
        output: cli.output,
        genes: cli.genes,
        mane_table: cli.mane_table,
        exclusions: Exclusions::new(cli.exclude_field),
        threads: cli.threads,
        command_line: std::env::args_os()
            .map(|arg| arg.to_string_lossy().into_owned())
            .collect::<Vec<_>>()
            .join(" "),
    };

    match convert(&config) {
        Ok(summary) => {
            info!(
                "wrote {} cells and {} variants to {}",
                summary.cells,
                summary.variants,
                config.output.display()
            );
            process::exit(0);
        }
        Err(e) => {
            error!("{e:#}");
            process::exit(1);
        }
    }
}

#[cfg(test)]
mod tests {
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

    #[test]
    fn help_text_is_plain_without_color() {
        let text = "Use `CSQ/SYMBOL`.\n  1  one copy\n";
        assert_eq!(
            style_help_text(text, false),
            "Use CSQ/SYMBOL.\n  1  one copy\n"
        );
    }
}
