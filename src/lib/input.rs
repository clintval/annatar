use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::num::NonZeroUsize;
use std::path::Path;

use anyhow::{Context, bail};
use noodles::bcf;
use noodles::bgzf;
use noodles::vcf::variant::RecordBuf;
use noodles::vcf::variant::record_buf::Samples;
use noodles::vcf::variant::record_buf::info::field::Value as InfoValue;
use noodles::vcf::variant::record_buf::info::field::value::Array as InfoArray;
use noodles::vcf::variant::record_buf::samples::sample::Value as SampleValue;
use noodles::vcf::variant::record_buf::samples::sample::value::Array as SampleArray;
use noodles::vcf::{self, Header};
use percent_encoding::percent_decode_str;

const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];
const BCF_MAGIC: [u8; 4] = *b"BCF\x02";

pub struct VariantSource {
    pub header: Header,
    pub header_lines: Vec<String>,
    pub samples: Vec<String>,
    pub records: Box<dyn Iterator<Item = io::Result<RecordBuf>>>,
}

// `Header` groups lines by key, which loses their order across keys, so the raw lines are kept.
struct HeaderTee {
    inner: Box<dyn BufRead + Send>,
    seen: Option<Vec<u8>>,
}

impl HeaderTee {
    fn new(inner: Box<dyn BufRead + Send>) -> Self {
        Self {
            inner,
            seen: Some(Vec::new()),
        }
    }

    fn take_lines(&mut self) -> Vec<String> {
        let seen = self.seen.take().unwrap_or_default();
        String::from_utf8_lossy(&seen)
            .lines()
            .filter(|line| line.starts_with("##"))
            .map(str::to_string)
            .collect()
    }
}

impl Read for HeaderTee {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.inner.read(buf)?;
        if let Some(seen) = &mut self.seen {
            seen.extend_from_slice(&buf[..n]);
        }
        Ok(n)
    }
}

impl BufRead for HeaderTee {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        self.inner.fill_buf()
    }

    fn consume(&mut self, amt: usize) {
        if let (Some(seen), Ok(buf)) = (&mut self.seen, self.inner.fill_buf()) {
            seen.extend_from_slice(&buf[..amt]);
        }
        self.inner.consume(amt);
    }
}

pub fn open_variants(path: &Path, threads: NonZeroUsize) -> anyhow::Result<VariantSource> {
    let raw: Box<dyn Read + Send> = if path == Path::new("-") {
        Box::new(io::stdin())
    } else {
        Box::new(File::open(path).with_context(|| format!("opening {}", path.display()))?)
    };
    let mut raw = BufReader::new(raw);

    let inner: Box<dyn BufRead + Send> = if raw.fill_buf()?.starts_with(&GZIP_MAGIC) {
        if !is_bgzf(raw.fill_buf()?) {
            let name = if path == Path::new("-") {
                "stdin".to_string()
            } else {
                path.display().to_string()
            };
            bail!("{name} is gzip but not BGZF; recompress it with bgzip");
        }
        let decoder = bgzf::io::MultithreadedReader::with_worker_count(threads, raw);
        Box::new(BufReader::new(decoder))
    } else {
        Box::new(raw)
    };

    let mut inner = HeaderTee::new(inner);
    if inner.fill_buf()?.starts_with(&BCF_MAGIC) {
        let mut reader = bcf::io::Reader::from(inner);
        let header = read_header(|| reader.read_header())?;
        let header_lines = reader.get_mut().take_lines();
        let samples = checked_samples(&header)?;
        let records = record_iter(
            header.clone(),
            move |record: &mut bcf::Record| reader.read_record(record),
            |record, header| {
                let name = record.reference_sequence_name(header.string_maps())?;
                Ok(format!("{name}:{}", position(record.variant_start())?))
            },
            true,
        );
        Ok(VariantSource {
            header,
            header_lines,
            samples,
            records: Box::new(records),
        })
    } else {
        let mut reader = vcf::io::Reader::new(inner);
        let header = read_header(|| reader.read_header())?;
        let header_lines = reader.get_mut().take_lines();
        let samples = checked_samples(&header)?;
        let records = record_iter(
            header.clone(),
            move |record: &mut vcf::Record| reader.read_record(record),
            |record, _| {
                let name = record.reference_sequence_name();
                Ok(format!("{name}:{}", position(record.variant_start())?))
            },
            false,
        );
        Ok(VariantSource {
            header,
            header_lines,
            samples,
            records: Box::new(records),
        })
    }
}

fn is_bgzf(head: &[u8]) -> bool {
    const FEXTRA: u8 = 0x04;
    if head.len() < 12 {
        return true;
    }
    if head[3] & FEXTRA == 0 {
        return false;
    }
    let extra_end = (12 + usize::from(u16::from_le_bytes([head[10], head[11]]))).min(head.len());
    let mut extra = &head[12..extra_end];
    while extra.len() >= 4 {
        if extra.starts_with(b"BC") {
            return true;
        }
        let len = usize::from(u16::from_le_bytes([extra[2], extra[3]]));
        extra = extra.get(4 + len..).unwrap_or_default();
    }
    false
}

fn position(start: Option<io::Result<noodles::core::Position>>) -> io::Result<usize> {
    match start {
        Some(result) => result.map(usize::from),
        None => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "missing position",
        )),
    }
}

fn record_iter<R>(
    header: Header,
    mut read: impl FnMut(&mut R) -> io::Result<usize>,
    locus_of: impl Fn(&R, &Header) -> io::Result<String>,
    needs_decoding: bool,
) -> impl Iterator<Item = io::Result<RecordBuf>>
where
    R: Default + vcf::variant::Record,
{
    let mut record = R::default();
    let mut ordinal = 0usize;
    std::iter::from_fn(move || {
        ordinal += 1;
        match read(&mut record) {
            Ok(0) => return None,
            Ok(_) => {}
            Err(e) => return Some(Err(with_context(&format!("record {ordinal}"), &e))),
        }
        let locus = match locus_of(&record, &header) {
            Ok(locus) => locus,
            Err(e) => return Some(Err(with_context(&format!("record {ordinal}"), &e))),
        };
        Some(
            RecordBuf::try_from_variant_record(&header, &record)
                .and_then(|mut buf| {
                    if needs_decoding {
                        normalize_bcf(&mut buf)?;
                    }
                    Ok(buf)
                })
                .map_err(|e| with_context(&locus, &e)),
        )
    })
}

fn with_context(locus: &str, e: &io::Error) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, format!("{locus}: {e}"))
}

fn read_header(read: impl FnOnce() -> io::Result<Header>) -> anyhow::Result<Header> {
    read().map_err(|e| {
        let message = e.to_string();
        match message.split_once("duplicate sample name: ") {
            Some((_, name)) => anyhow::anyhow!("sample {name} appears more than once"),
            None => anyhow::anyhow!("reading the VCF header: {message}"),
        }
    })
}

fn checked_samples(header: &Header) -> anyhow::Result<Vec<String>> {
    if header.sample_names().is_empty() {
        bail!("the VCF has no samples");
    }
    Ok(header.sample_names().iter().cloned().collect())
}

fn decode(s: &str) -> io::Result<String> {
    if !s.contains('%') {
        return Ok(s.to_string());
    }
    percent_decode_str(s)
        .decode_utf8()
        .map(std::borrow::Cow::into_owned)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn is_missing(array: &SampleArray) -> bool {
    match array {
        SampleArray::Integer(v) => matches!(v.as_slice(), [None]),
        SampleArray::Float(v) => matches!(v.as_slice(), [None]),
        SampleArray::Character(v) => matches!(v.as_slice(), [None]),
        SampleArray::String(v) => matches!(v.as_slice(), [None]),
    }
}

fn is_missing_info(value: &InfoValue) -> bool {
    match value {
        InfoValue::Character(c) => *c == '.',
        InfoValue::String(s) => s == ".",
        InfoValue::Array(InfoArray::Character(v)) => matches!(v.as_slice(), [None]),
        InfoValue::Array(InfoArray::String(v)) => matches!(v.as_slice(), [None]),
        _ => false,
    }
}

fn normalize_bcf(record: &mut RecordBuf) -> io::Result<()> {
    let info = record.info_mut();
    for i in 0..info.keys().count() {
        if let Some((_, value)) = info.get_index_mut(i) {
            if value.as_ref().is_some_and(is_missing_info) {
                *value = None;
            } else if let Some(InfoValue::String(s)) = value {
                *s = decode(s)?;
            }
        }
    }

    let keys = record.samples().keys().clone();
    let values = record
        .samples()
        .values()
        .map(|sample| {
            sample
                .values()
                .iter()
                .map(|value| match value {
                    Some(SampleValue::String(s)) => Ok(Some(SampleValue::String(decode(s)?))),
                    Some(SampleValue::Array(a)) if is_missing(a) => Ok(None),
                    other => Ok(other.clone()),
                })
                .collect::<io::Result<Vec<_>>>()
        })
        .collect::<io::Result<Vec<_>>>()?;
    *record.samples_mut() = Samples::new(keys, values);
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::path::PathBuf;

    use noodles::vcf::variant::record_buf::info::field::Value;
    use noodles::vcf::variant::record_buf::info::field::value::Array;
    use rstest::rstest;

    use super::*;

    fn data(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data")
            .join(name)
    }

    fn open(path: &Path) -> anyhow::Result<VariantSource> {
        open_variants(path, NonZeroUsize::new(1).unwrap())
    }

    fn collect(name: &str) -> Vec<RecordBuf> {
        open(&data(name))
            .unwrap()
            .records
            .collect::<io::Result<_>>()
            .unwrap()
    }

    fn csq_strings(record: &RecordBuf) -> Vec<String> {
        match record.info().get("CSQ").flatten() {
            Some(Value::Array(Array::String(v))) => v.iter().flatten().cloned().collect(),
            other => panic!("unexpected CSQ {other:?}"),
        }
    }

    #[rstest]
    fn keeps_the_raw_header_lines(#[values("basic.vcf", "basic.vcf.gz", "basic.bcf")] name: &str) {
        let lines = open(&data(name)).unwrap().header_lines;
        assert!(lines.iter().all(|line| line.starts_with("##")));
        let vep = |prefix: &str| lines.iter().position(|line| line.starts_with(prefix));
        assert!(vep("##VEP=").unwrap() < vep("##VEP-command-line=").unwrap());
    }

    #[rstest]
    fn opens_every_format(#[values("basic.vcf", "basic.vcf.gz", "basic.bcf")] name: &str) {
        let source = open(&data(name)).unwrap();
        assert_eq!(source.samples, ["c1", "c2", "c3"]);
        let records = collect(name);
        assert_eq!(records.len(), 4);
        let expected = collect("basic.vcf");
        for (a, b) in records.iter().zip(&expected) {
            assert_eq!(a.reference_sequence_name(), b.reference_sequence_name());
            assert_eq!(a.variant_start(), b.variant_start());
            assert_eq!(a.reference_bases(), b.reference_bases());
            assert_eq!(a.alternate_bases(), b.alternate_bases());
            assert_eq!(a.quality_score(), b.quality_score());
            assert_eq!(a.filters(), b.filters());
            assert_eq!(a.info(), b.info());
            assert_eq!(a.samples(), b.samples());
        }
    }

    #[rstest]
    fn decodes_percent_encoding(#[values("basic.vcf", "basic.vcf.gz", "basic.bcf")] name: &str) {
        let records = collect(name);
        let csq = csq_strings(&records[1]);
        assert!(
            csq.iter().any(|s| s.contains("Li-Fraumeni,syndrome")),
            "{csq:?}"
        );
        assert!(!csq.iter().any(|s| s.contains("%2C")));
    }

    #[rstest]
    fn decodes_exactly_once(#[values("basic.vcf", "basic.vcf.gz", "basic.bcf")] name: &str) {
        let records = collect(name);
        let csq = csq_strings(&records[2]);
        assert!(csq.iter().any(|s| s.contains("Odd%25name")), "{csq:?}");
    }

    #[test]
    fn missing_arrays_become_none() {
        let mut record = RecordBuf::default();
        let keys = ["AD".to_string()].into_iter().collect();
        let values = vec![
            vec![Some(SampleValue::Array(SampleArray::Integer(vec![None])))],
            vec![Some(SampleValue::Array(SampleArray::Integer(vec![
                Some(1),
                None,
            ])))],
        ];
        *record.samples_mut() = Samples::new(keys, values);
        normalize_bcf(&mut record).unwrap();
        let got: Vec<_> = record
            .samples()
            .values()
            .map(|s| s.values()[0].clone())
            .collect();
        assert_eq!(got[0], None);
        assert_eq!(
            got[1],
            Some(SampleValue::Array(SampleArray::Integer(vec![
                Some(1),
                None
            ])))
        );
    }

    #[test]
    fn decodes_scalar_strings_once() {
        let mut record = RecordBuf::default();
        record
            .info_mut()
            .insert("S".into(), Some(Value::String("a%2525b".into())));
        normalize_bcf(&mut record).unwrap();
        assert_eq!(
            record.info().get("S").flatten(),
            Some(&Value::String("a%25b".into()))
        );
    }

    #[test]
    fn missing_scalar_info_becomes_none() {
        let mut record = RecordBuf::default();
        record
            .info_mut()
            .insert("S".into(), Some(Value::String(".".into())));
        record
            .info_mut()
            .insert("C".into(), Some(Value::Character('.')));
        normalize_bcf(&mut record).unwrap();
        assert_eq!(record.info().get("S"), Some(None));
        assert_eq!(record.info().get("C"), Some(None));
    }

    #[test]
    fn a_lone_missing_text_array_becomes_none() {
        let partial = Value::Array(Array::String(vec![Some("x".into()), None]));
        let mut record = RecordBuf::default();
        let info = record.info_mut();
        info.insert("S".into(), Some(Value::Array(Array::String(vec![None]))));
        info.insert("C".into(), Some(Value::Array(Array::Character(vec![None]))));
        info.insert("P".into(), Some(partial.clone()));
        normalize_bcf(&mut record).unwrap();
        assert_eq!(record.info().get("S"), Some(None));
        assert_eq!(record.info().get("C"), Some(None));
        assert_eq!(record.info().get("P"), Some(Some(&partial)));
    }

    #[test]
    fn names_the_record_on_a_malformed_line() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("bad.vcf");
        let text = "##fileformat=VCFv4.3\n##contig=<ID=chr1>\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tc1\n\
            chr1\tabc\t.\tA\tG\t.\t.\t.\tGT\t0/1\n";
        fs::write(&path, text).unwrap();
        let err = open(&path).unwrap().records.next().unwrap().unwrap_err();
        assert!(err.to_string().contains("record 1"), "{err}");
    }

    #[test]
    fn rejects_sites_only() {
        let err = open(&data("sites-only.vcf")).err().unwrap();
        assert!(err.to_string().contains("the VCF has no samples"), "{err}");
    }

    #[test]
    fn rejects_duplicate_samples() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("dup.vcf");
        let text = fs::read_to_string(data("basic.vcf"))
            .unwrap()
            .replace("\tc2\tc3", "\tc1\tc3");
        fs::write(&path, text).unwrap();
        let err = open(&path).err().unwrap();
        assert!(
            err.to_string().contains("sample c1 appears more than once"),
            "{err}"
        );
    }

    #[test]
    fn names_the_record_on_a_bad_value() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("bad.vcf");
        let text = "##fileformat=VCFv4.3\n##contig=<ID=chr1>\n\
            ##INFO=<ID=N,Number=1,Type=Integer,Description=\"n\">\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tc1\n\
            chr1\t42\t.\tA\tG\t.\t.\tN=abc\tGT\t0/1\n";
        fs::write(&path, text).unwrap();
        let err = open(&path).unwrap().records.next().unwrap().unwrap_err();
        assert!(err.to_string().contains("chr1:42"), "{err}");
    }
}
