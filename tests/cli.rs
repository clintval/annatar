use std::fs;
use std::path::{Path, PathBuf};

use anndata::{AnnData, AnnDataOp, Backend};
use anndata_hdf5::H5;
use assert_cmd::Command;
use predicates::prelude::PredicateBooleanExt;
use predicates::str::contains;
use rstest::rstest;
use tempfile::TempDir;

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

fn annatar() -> Command {
    let mut cmd = Command::cargo_bin("annatar").unwrap();
    cmd.env("NO_COLOR", "1").env_remove("RUST_LOG");
    cmd
}

fn output_in(dir: &TempDir) -> PathBuf {
    dir.path().join("out.h5ad")
}

fn open(path: &Path) -> AnnData<H5> {
    AnnData::<H5>::open(H5::open(path).unwrap()).unwrap()
}

#[test]
fn converts_basic() {
    let dir = TempDir::new().unwrap();
    let output = output_in(&dir);
    annatar()
        .arg("--output")
        .arg(&output)
        .arg(data("basic.vcf"))
        .assert()
        .success();
    assert!(output.exists());
}

#[test]
fn reads_stdin() {
    let dir = TempDir::new().unwrap();
    let output = output_in(&dir);
    annatar()
        .arg("-o")
        .arg(&output)
        .arg("-")
        .write_stdin(fs::read(data("basic.vcf")).unwrap())
        .assert()
        .success();
    let adata = open(&output);
    assert_eq!((adata.n_obs(), adata.n_vars()), (3, 5));
}

#[test]
fn stdin_is_the_default_input() {
    let dir = TempDir::new().unwrap();
    let output = output_in(&dir);
    annatar()
        .arg("-o")
        .arg(&output)
        .write_stdin(fs::read(data("basic.vcf")).unwrap())
        .assert()
        .success();
    assert!(output.exists());
}

#[test]
fn threads_read_a_bgzipped_input() {
    let dir = TempDir::new().unwrap();
    let output = output_in(&dir);
    annatar()
        .arg("-o")
        .arg(&output)
        .args(["-@", "2"])
        .arg(data("basic.vcf.gz"))
        .assert()
        .success();
    assert_eq!(open(&output).n_vars(), 5);
}

#[test]
fn plain_gzip_input_is_refused_with_advice() {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("plain.vcf.gz");
    let mut encoder = flate2::write::GzEncoder::new(
        fs::File::create(&input).unwrap(),
        flate2::Compression::default(),
    );
    std::io::Write::write_all(&mut encoder, &fs::read(data("basic.vcf")).unwrap()).unwrap();
    encoder.finish().unwrap();
    let output = output_in(&dir);
    annatar()
        .arg("--output")
        .arg(&output)
        .arg(&input)
        .assert()
        .failure()
        .stderr(contains(
            "plain.vcf.gz is gzip but not BGZF; recompress it with bgzip",
        ));
    assert!(!output.exists());
}

#[test]
fn no_records_writes_empty_matrix() {
    let dir = TempDir::new().unwrap();
    let output = output_in(&dir);
    annatar()
        .arg("-o")
        .arg(&output)
        .arg(data("no-records.vcf"))
        .assert()
        .success();
    let adata = open(&output);
    assert_eq!((adata.n_obs(), adata.n_vars()), (3, 0));
    assert_eq!(adata.obs_names().into_vec(), ["c1", "c2", "c3"]);
}

#[test]
fn symbolic_alts_are_skipped() {
    let dir = TempDir::new().unwrap();
    let output = output_in(&dir);
    annatar()
        .arg("-o")
        .arg(&output)
        .arg(data("symbolic.vcf"))
        .assert()
        .success()
        .stderr(contains("symbolic ALT alleles skipped: 2").count(1));

    let adata = open(&output);
    assert_eq!(
        adata.var_names().into_vec(),
        [
            "chr1:100:AC:AT",
            "chr1:200:A:G",
            "chr1:300:A:G",
            "chr1:300:A:*"
        ]
    );
    let var = adata.read_var().unwrap();
    let symbol = var.column("SYMBOL").unwrap();
    assert_eq!(symbol.null_count(), 1);
    assert!(symbol.get(3).unwrap().is_null());
    assert!(var.column("Consequence").unwrap().get(3).unwrap().is_null());
    assert!(var.column("pick_tier").unwrap().get(3).unwrap().is_null());
}

#[rstest]
#[case("no-ad.vcf", "", "the VCF declares no FORMAT/AD")]
#[case(
    "ad-number-dot.vcf",
    "",
    "FORMAT/AD must be declared Number=R, Type=Integer"
)]
#[case(
    "ad-length.vcf",
    "chr1:100",
    "FORMAT/AD has 3 values but this record needs 2"
)]
#[case(
    "ad-length-uncovered.vcf",
    "chr1:100",
    "FORMAT/AD has 3 values but this record needs 2"
)]
#[case(
    "ad-length-counts-symbolic.vcf",
    "chr1:100",
    "FORMAT/AD has 2 values but this record needs 3"
)]
#[case(
    "csq-subfields.vcf",
    "chr1:100",
    "a CSQ entry has 3 subfields where the header declares 7"
)]
#[case("no-csq-entry.vcf", "chr1:100", "no CSQ entry for ALT T")]
#[case("no-csq-entry-symbolic.vcf", "chr1:100", "no CSQ entry for ALT <DEL>")]
#[case(
    "no-csq-entry-breakend.vcf",
    "chr1:100",
    "no CSQ entry for ALT A[chr1:500["
)]
#[case(
    "layer-collides-with-ad.vcf",
    "",
    "FORMAT/AD and FORMAT/AD_ref both make the layer AD_ref"
)]
#[case(
    "layer-collides-number-r.vcf",
    "",
    "FORMAT/X and FORMAT/X_ref both make the layer X_ref"
)]
#[case(
    "layer-collides-fixed.vcf",
    "",
    "FORMAT/X and FORMAT/X_1 both make the layer X_1"
)]
#[case(
    "var-column-collides-info.vcf",
    "",
    "two fields make the var column X_ref"
)]
#[case(
    "var-column-collides-csq.vcf",
    "",
    "two fields make the var column pos"
)]
#[case(
    "csq-no-format.vcf",
    "",
    "INFO/CSQ has no \"Format: \" list in its Description"
)]
#[case(
    "gt-allele-range.vcf",
    "chr1:100",
    "FORMAT/GT calls allele 3 but the last allele here is 1"
)]
#[case(
    "negative-ad-uncovered.vcf",
    "chr1:100",
    "FORMAT/AD holds a negative count, -5"
)]
#[case(
    "negative-ad-skipped-alt.vcf",
    "chr1:100",
    "FORMAT/AD holds a negative count, -1"
)]
#[case(
    "unranked-consequence.vcf",
    "chr1:100",
    "does not rank the Consequence term 'invented_variant'"
)]
#[case("bad-info-value.vcf", "chr1:100", "invalid digit")]
#[case("bad-format-value.vcf", "chr1:100", "invalid digit")]
#[case("duplicate-sample.vcf", "", "sample c1 appears more than once")]
#[case(
    "duplicate-variant.vcf",
    "chr1:100",
    "variant chr1:100:A:G appears more than once"
)]
#[case(
    "info-number-a-length.vcf",
    "chr1:100",
    "INFO/AC has 2 values but this record needs 1"
)]
#[case(
    "info-number-r-length.vcf",
    "chr1:100",
    "INFO/RD has 1 values but this record needs 2"
)]
#[case(
    "info-fixed-length.vcf",
    "chr1:100",
    "INFO/QUARTET has 3 values but this record needs 4"
)]
#[case("negative-ad.vcf", "chr1:100", "FORMAT/AD holds a negative count, -1")]
#[case(
    "format-number-a-length.vcf",
    "chr1:100",
    "FORMAT/AF has 2 values but this record needs 1"
)]
fn an_error_input_fails_naming_the_problem_and_writes_nothing(
    #[case] name: &str,
    #[case] locus: &str,
    #[case] message: &str,
) {
    let dir = TempDir::new().unwrap();
    let output = output_in(&dir);
    annatar()
        .arg("-o")
        .arg(&output)
        .arg(data("errors").join(name))
        .assert()
        .failure()
        .stderr(contains(locus).and(contains(message)));
    assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 0);
}

#[rstest]
#[case("basic.vcf", "FORMAT/AD", "FORMAT/AD is required")]
#[case("basic.vcf", "FORMAT/GT", "FORMAT/GT cannot be excluded")]
#[case("basic.vcf", "INFO/CSQ", "INFO/CSQ cannot be excluded")]
#[case(
    "basic.vcf",
    "CSQ/symbol",
    "--exclude-field CSQ/symbol names nothing in this VCF's header"
)]
#[case(
    "basic.vcf",
    "INFO/SYMBOL",
    "--exclude-field INFO/SYMBOL names nothing in this VCF's header"
)]
#[case(
    "basic.vcf",
    "FORMAT/dp",
    "--exclude-field FORMAT/dp names nothing in this VCF's header"
)]
#[case(
    "no-csq.vcf",
    "CSQ/SYMBOL",
    "--exclude-field CSQ/SYMBOL names nothing in this VCF's header"
)]
fn an_exclusion_that_cannot_be_honored_is_an_error(
    #[case] name: &str,
    #[case] field: &str,
    #[case] message: &str,
) {
    let dir = TempDir::new().unwrap();
    annatar()
        .arg("-o")
        .arg(output_in(&dir))
        .args(["--exclude-field", field])
        .arg(data(name))
        .assert()
        .failure()
        .stderr(contains(message));
    assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 0);
}

#[test]
fn a_malformed_exclude_field_is_an_error() {
    let dir = TempDir::new().unwrap();
    annatar()
        .arg("-o")
        .arg(output_in(&dir))
        .args(["--exclude-field", "SYMBOL"])
        .arg(data("basic.vcf"))
        .assert()
        .failure()
        .stderr(contains("expected INFO/ID, FORMAT/ID or CSQ/subfield"));
}

#[rstest]
#[case("--genes", "no-such-genes.txt")]
#[case("--mane-table", "no-such-mane.txt")]
fn a_missing_lookup_file_is_an_error_naming_the_path(#[case] flag: &str, #[case] file: &str) {
    let dir = TempDir::new().unwrap();
    let missing = dir.path().join(file);
    annatar()
        .arg("-o")
        .arg(output_in(&dir))
        .arg(flag)
        .arg(&missing)
        .arg(data("basic.vcf"))
        .assert()
        .failure()
        .stderr(contains(missing.display().to_string()));
}

#[test]
fn an_unwritable_output_is_an_error_naming_the_path() {
    let dir = TempDir::new().unwrap();
    let output = dir.path().join("no-such-directory").join("out.h5ad");
    annatar()
        .arg("-o")
        .arg(&output)
        .arg(data("basic.vcf"))
        .assert()
        .failure()
        .stderr(contains(output.display().to_string()));
}

#[test]
fn a_failed_run_keeps_an_earlier_output() {
    let dir = TempDir::new().unwrap();
    let output = output_in(&dir);
    fs::write(&output, "earlier").unwrap();
    annatar()
        .arg("-o")
        .arg(&output)
        .arg(data("errors/ad-length.vcf"))
        .assert()
        .failure();
    assert_eq!(fs::read_to_string(&output).unwrap(), "earlier");
}

#[test]
fn a_successful_run_replaces_an_earlier_output() {
    let dir = TempDir::new().unwrap();
    let output = output_in(&dir);
    fs::write(&output, "earlier").unwrap();
    annatar()
        .arg("-o")
        .arg(&output)
        .arg(data("basic.vcf"))
        .assert()
        .success();
    assert_eq!(open(&output).n_vars(), 5);
    assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 1);
}

#[test]
fn a_failed_write_keeps_an_earlier_output() {
    let dir = TempDir::new().unwrap();
    let output = output_in(&dir);
    fs::write(&output, "earlier").unwrap();
    fs::create_dir(dir.path().join("out.h5ad.partial")).unwrap();
    annatar()
        .arg("-o")
        .arg(&output)
        .arg(data("basic.vcf"))
        .assert()
        .failure()
        .stderr(contains(output.display().to_string()));
    assert_eq!(fs::read_to_string(&output).unwrap(), "earlier");
}

#[test]
fn a_failed_rename_leaves_no_partial_file() {
    let dir = TempDir::new().unwrap();
    let output = output_in(&dir);
    fs::create_dir(&output).unwrap();
    annatar()
        .arg("-o")
        .arg(&output)
        .arg(data("basic.vcf"))
        .assert()
        .failure()
        .stderr(contains(output.display().to_string()));
    assert!(!dir.path().join("out.h5ad.partial").exists());
}

#[test]
fn output_names_differing_by_suffix_do_not_share_a_partial_file() {
    let dir = TempDir::new().unwrap();
    fs::create_dir(dir.path().join("cells.h5ad.partial")).unwrap();
    let output = dir.path().join("cells.1");
    annatar()
        .arg("-o")
        .arg(&output)
        .arg(data("basic.vcf"))
        .assert()
        .success();
    assert_eq!(open(&output).n_vars(), 5);
}

#[test]
fn standard_output_is_refused() {
    let dir = TempDir::new().unwrap();
    annatar()
        .current_dir(dir.path())
        .args(["-o", "-"])
        .arg(data("basic.vcf"))
        .assert()
        .failure()
        .stderr(contains("HDF5 cannot stream"));
    assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 0);
}

#[cfg(unix)]
#[test]
fn a_non_utf8_path_is_an_ordinary_error() {
    use std::ffi::OsStr;
    use std::os::unix::ffi::OsStrExt;

    let dir = TempDir::new().unwrap();
    annatar()
        .arg("-o")
        .arg(output_in(&dir))
        .arg(OsStr::from_bytes(b"in\xff.vcf"))
        .assert()
        .code(1)
        .stderr(contains("opening"));
}

#[test]
fn records_with_no_alt_warn_once() {
    let dir = TempDir::new().unwrap();
    let output = output_in(&dir);
    annatar()
        .arg("-o")
        .arg(&output)
        .arg(data("no-alt.vcf"))
        .assert()
        .success()
        .stderr(contains("records with no ALT skipped: 1").count(1));
    assert_eq!(open(&output).n_vars(), 1);
}

#[test]
fn a_mane_table_is_ignored_with_a_warning_when_the_header_declares_mane_fields() {
    let dir = TempDir::new().unwrap();
    annatar()
        .arg("-o")
        .arg(output_in(&dir))
        .arg("--mane-table")
        .arg(data("mane-summary.txt"))
        .arg(data("basic.vcf"))
        .assert()
        .success()
        .stderr(contains("--mane-table is ignored").count(1));
}

#[test]
fn no_gt_warns_once() {
    let dir = TempDir::new().unwrap();
    annatar()
        .arg("-o")
        .arg(output_in(&dir))
        .arg(data("no-gt.vcf"))
        .assert()
        .success()
        .stderr(contains("no FORMAT/GT").count(1));
}

#[test]
fn missing_tie_fields_warn_once_with_csq() {
    let dir = TempDir::new().unwrap();
    annatar()
        .arg("-o")
        .arg(output_in(&dir))
        .arg(data("multiallelic.vcf"))
        .assert()
        .success()
        .stderr(contains("MANE_SELECT").count(1));
}

#[test]
fn missing_tie_fields_stay_quiet_without_csq() {
    let dir = TempDir::new().unwrap();
    annatar()
        .arg("-o")
        .arg(output_in(&dir))
        .arg(data("no-csq.vcf"))
        .assert()
        .success()
        .stderr(contains("MANE_SELECT").not());
}

#[test]
fn missing_values_warn_once_per_field() {
    let dir = TempDir::new().unwrap();
    annatar()
        .arg("-o")
        .arg(output_in(&dir))
        .arg(data("missing-ref-count.vcf"))
        .assert()
        .success()
        .stderr(contains("AD_ref: missing in 1 of the covered cells").count(1));
}

#[test]
fn a_clean_run_has_no_warnings() {
    let dir = TempDir::new().unwrap();
    annatar()
        .arg("-o")
        .arg(output_in(&dir))
        .arg(data("basic.vcf"))
        .assert()
        .success()
        .stderr(contains("WARN").not());
}

#[test]
fn the_help_lists_every_flag() {
    annatar().arg("--help").assert().success().stdout(
        contains("--output")
            .and(contains("--genes"))
            .and(contains("--mane-table"))
            .and(contains("--exclude-field"))
            .and(contains("--threads"))
            .and(contains("MIT License"))
            .and(contains("`").not()),
    );
}
