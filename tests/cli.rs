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
        .stderr(contains("*, <*> and <NON_REF> ALT alleles skipped: 3").count(1));

    let adata = open(&output);
    assert_eq!(
        adata.var_names().into_vec(),
        ["chr1:100:AC:AT", "chr1:200:A:G", "chr1:300:A:G"]
    );
}

#[test]
fn a_structural_alt_with_no_entry_keeps_its_row_and_warns_once() {
    let dir = TempDir::new().unwrap();
    let output = output_in(&dir);
    annatar()
        .arg("-o")
        .arg(&output)
        .arg(data("structural-no-entry.vcf"))
        .assert()
        .success()
        .stderr(contains("symbolic or breakend ALT alleles with no CSQ entry: 2").count(1));
    assert_eq!(open(&output).n_vars(), 4);
}

#[rstest]
#[case("no-csq-entry.vcf")]
#[case("no-csq-entry-allele-number.vcf")]
#[case("split-after-vep.vcf")]
#[case("split-after-vep-allele-number.vcf")]
#[case("pick-then-split.vcf")]
#[case("single-alt-mismatch.vcf")]
fn allele_number_is_not_suggested_where_it_would_not_help(#[case] name: &str) {
    let dir = TempDir::new().unwrap();
    annatar()
        .arg("-o")
        .arg(output_in(&dir))
        .arg(data("errors").join(name))
        .assert()
        .failure()
        .stderr(contains("--allele_number").not());
}

#[test]
fn a_minimal_breakend_with_no_number_keeps_an_empty_row_and_warns_once() {
    let dir = TempDir::new().unwrap();
    let output = output_in(&dir);
    annatar()
        .arg("-o")
        .arg(&output)
        .arg(data("minimal-no-number.vcf"))
        .assert()
        .success()
        .stderr(contains("symbolic or breakend ALT alleles with no CSQ entry: 1").count(1));
    assert_eq!(open(&output).n_vars(), 3);
}

#[rstest]
#[case("symbolic-csq.vcf")]
#[case("breakend-csq.vcf")]
#[case("breakend-reference-side.vcf")]
#[case("allele-number.vcf")]
#[case("allele-spelling.vcf")]
#[case("split-after-vep-snv.vcf")]
#[case("split-after-vep-breakend.vcf")]
fn alleles_vep_spells_its_own_way_match_without_warning(#[case] name: &str) {
    let dir = TempDir::new().unwrap();
    annatar()
        .arg("-o")
        .arg(output_in(&dir))
        .arg(data(name))
        .assert()
        .success()
        .stderr(contains("no CSQ entry").not());
}

#[rstest]
#[case(
    "no-ad.vcf",
    "",
    "the VCF declares no FORMAT/AD; with bcftools mpileup, add --annotate FORMAT/AD"
)]
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
#[case(
    "no-csq-entry.vcf",
    "chr1:100",
    "no CSQ entry for ALT T among the Allele values G, C; annotate after splitting"
)]
#[case(
    "no-csq-entry-allele-number.vcf",
    "chr1:100",
    "no CSQ entry for ALT T among the Allele values G; annotate after splitting"
)]
#[case(
    "split-after-vep.vcf",
    "chr1:100",
    "no CSQ entry for ALT AT, spelled T by VEP, among the Allele values G, AT; \
     annotate after splitting"
)]
#[case(
    "split-after-vep-allele-number.vcf",
    "chr1:100",
    "no CSQ entry for ALT AT, spelled T by VEP, among the Allele values G, AT; \
     annotate after splitting"
)]
#[case(
    "pick-then-split.vcf",
    "chr1:100",
    "no CSQ entry for ALT T among the Allele values A; annotate after splitting"
)]
#[case(
    "single-alt-mismatch.vcf",
    "chr1:100",
    "no CSQ entry for ALT G among the Allele values T; annotate after splitting"
)]
#[case(
    "minimal-multiallelic.vcf",
    "chr1:100",
    "no CSQ entry for ALT ATT, spelled TT by VEP, among the Allele values -, T; \
     annotate after splitting, or with VEP's --allele_number for multiallelic --minimal output"
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
#[case("bad-info-value.vcf", "chr1:100", "invalid digit")]
#[case("bad-format-value.vcf", "chr1:100", "invalid digit")]
#[case("duplicate-sample.vcf", "", "sample c1 appears more than once")]
#[case(
    "duplicate-variant.vcf",
    "chr1:100",
    "variant chr1:100:A:G appears more than once; remove repeats with bcftools norm -d exact"
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
#[case(
    "basic.vcf",
    "INFO/CSQ/symbol",
    "--exclude-field INFO/CSQ/symbol names nothing in this VCF's header"
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
    "INFO/SYMBOL",
    "--exclude-field INFO/SYMBOL names nothing in this VCF's header"
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

const NO_CSQ_NOTE: &str = "no INFO/CSQ to read, so the variants have no annotation";

#[test]
fn a_header_with_no_csq_gets_one_info_line_and_no_warning() {
    let dir = TempDir::new().unwrap();
    annatar()
        .arg("-o")
        .arg(output_in(&dir))
        .arg(data("no-csq.vcf"))
        .assert()
        .success()
        .stderr(contains(NO_CSQ_NOTE).count(1).and(contains("WARN").not()));
}

#[test]
fn a_header_with_csq_gets_no_note() {
    let dir = TempDir::new().unwrap();
    annatar()
        .arg("-o")
        .arg(output_in(&dir))
        .arg(data("basic.vcf"))
        .assert()
        .success()
        .stderr(contains(NO_CSQ_NOTE).not());
}

#[rstest]
#[case("--genes", "panel-genes.txt")]
#[case("--mane-table", "mane-summary.txt")]
fn a_lookup_file_warns_when_there_is_no_csq_to_act_on(#[case] flag: &str, #[case] file: &str) {
    let dir = TempDir::new().unwrap();
    annatar()
        .arg("-o")
        .arg(output_in(&dir))
        .arg(flag)
        .arg(data(file))
        .arg(data("no-csq.vcf"))
        .assert()
        .success()
        .stderr(
            contains(format!(
                "{flag} does nothing: there is no INFO/CSQ to act on"
            ))
            .count(1),
        );
}

#[rstest]
#[case("no-csq.vcf", "INFO/CSQ/SYMBOL")]
#[case("no-csq.vcf", "INFO/CSQ")]
fn an_exclusion_of_absent_annotation_warns(#[case] name: &str, #[case] field: &str) {
    let dir = TempDir::new().unwrap();
    annatar()
        .arg("-o")
        .arg(output_in(&dir))
        .args(["--exclude-field", field])
        .arg(data(name))
        .assert()
        .success()
        .stderr(
            contains(format!(
                "--exclude-field {field} does nothing: there is no INFO/CSQ to act on"
            ))
            .count(1),
        );
}

#[test]
fn excluding_info_csq_acts_as_a_header_with_no_csq() {
    let dir = TempDir::new().unwrap();
    annatar()
        .arg("-o")
        .arg(output_in(&dir))
        .args(["--exclude-field", "INFO/CSQ"])
        .args(["--exclude-field", "INFO/CSQ/no_such_subfield"])
        .arg("--genes")
        .arg(data("panel-genes.txt"))
        .arg(data("basic.vcf"))
        .assert()
        .success()
        .stderr(
            contains(NO_CSQ_NOTE)
                .count(1)
                .and(contains("--exclude-field INFO/CSQ does nothing").not())
                .and(contains("--exclude-field INFO/CSQ/no_such_subfield does nothing").count(1))
                .and(contains("--genes does nothing").count(1)),
        );
}

#[test]
fn excluding_info_csq_allows_a_csq_header_with_no_format_list() {
    let dir = TempDir::new().unwrap();
    annatar()
        .arg("-o")
        .arg(output_in(&dir))
        .args(["--exclude-field", "INFO/CSQ"])
        .arg(data("errors/csq-no-format.vcf"))
        .assert()
        .success();
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
        .stderr(contains("expected INFO/ID, FORMAT/ID or INFO/CSQ/subfield"));
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
fn a_value_failing_its_csq_type_warns_once_per_column() {
    let dir = TempDir::new().unwrap();
    annatar()
        .arg("-o")
        .arg(output_in(&dir))
        .arg(data("csq-types-fallback.vcf"))
        .assert()
        .success()
        .stderr(
            contains("INFO/CSQ/STRAND: plus is not an integer, so the column is text in this file")
                .count(1)
                .and(contains("INFO/CSQ/TSL: 1&2 is not an integer").count(1))
                .and(contains("INFO/CSQ/REVEL_score: nan is not a number").count(1))
                .and(contains("INFO/CSQ/gnomADe_AF: inf is not a number").count(1))
                .and(contains("INFO/CSQ/CADD_PHRED: 1.5abc is not a number").count(1)),
        );
}

#[test]
fn an_unranked_consequence_warns_once_per_term() {
    let dir = TempDir::new().unwrap();
    annatar()
        .arg("-o")
        .arg(output_in(&dir))
        .arg(data("unranked-consequence.vcf"))
        .assert()
        .success()
        .stderr(
            contains("VEP 116 does not rank the Consequence term invented_variant")
                .count(1)
                .and(contains("Consequence term other_invented_variant").count(1)),
        );
}

#[test]
fn another_vep_major_warns_once_and_runs() {
    let dir = TempDir::new().unwrap();
    let output = output_in(&dir);
    annatar()
        .arg("-o")
        .arg(&output)
        .arg(data("vep-115.vcf"))
        .assert()
        .success()
        .stderr(contains("annotated by VEP 115").count(1));
    assert_eq!(open(&output).n_vars(), 5);
}

#[test]
fn a_header_with_no_vep_line_gets_no_version_warning() {
    let dir = TempDir::new().unwrap();
    annatar()
        .arg("-o")
        .arg(output_in(&dir))
        .arg(data("multiallelic.vcf"))
        .assert()
        .success()
        .stderr(contains("annotated by VEP").not());
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

const ESCAPE: &str = "\u{1b}";

fn annatar_with_no_color(value: Option<&str>) -> Command {
    let mut cmd = Command::cargo_bin("annatar").unwrap();
    cmd.env_remove("RUST_LOG").env_remove("NO_COLOR");
    if let Some(value) = value {
        cmd.env("NO_COLOR", value);
    }
    cmd
}

#[rstest]
#[case(None)]
#[case(Some(""))]
#[case(Some("1"))]
fn piped_help_has_no_escapes(#[case] no_color: Option<&str>) {
    annatar_with_no_color(no_color)
        .arg("--help")
        .assert()
        .success()
        .stdout(contains("--output").and(contains(ESCAPE).not()));
}

#[rstest]
#[case(None)]
#[case(Some(""))]
#[case(Some("1"))]
fn a_piped_usage_error_has_no_escapes(#[case] no_color: Option<&str>) {
    annatar_with_no_color(no_color)
        .arg("--bogus")
        .assert()
        .failure()
        .stderr(contains("--bogus").and(contains(ESCAPE).not()));
}

#[rstest]
#[case(None)]
#[case(Some(""))]
#[case(Some("1"))]
fn piped_help_for_no_arguments_has_no_escapes(#[case] no_color: Option<&str>) {
    annatar_with_no_color(no_color)
        .assert()
        .failure()
        .stderr(contains("--output").and(contains(ESCAPE).not()));
}

#[rstest]
#[case(None)]
#[case(Some(""))]
#[case(Some("1"))]
fn a_failing_run_has_no_escapes_in_a_pipe(#[case] no_color: Option<&str>) {
    let dir = TempDir::new().unwrap();
    annatar_with_no_color(no_color)
        .arg("-o")
        .arg(output_in(&dir))
        .arg(data("errors/ad-length.vcf"))
        .assert()
        .failure()
        .stderr(contains("ERROR").and(contains(ESCAPE).not()));
}

#[test]
fn forcing_color_in_the_environment_does_not_color_a_pipe() {
    annatar_with_no_color(None)
        .env("CLICOLOR_FORCE", "1")
        .arg("--bogus")
        .assert()
        .failure()
        .stderr(contains(ESCAPE).not());
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
            .and(contains("VEP-annotated multi-sample").not())
            .and(contains("`").not()),
    );
}
