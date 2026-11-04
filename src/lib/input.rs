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
    pub samples: Vec<String>,
    pub records: Box<dyn Iterator<Item = io::Result<RecordBuf>>>,
}

pub fn open_variants(path: &Path, threads: NonZeroUsize) -> anyhow::Result<VariantSource> {
    let raw: Box<dyn Read + Send> = if path == Path::new("-") {
        Box::new(io::stdin())
    } else {
        Box::new(File::open(path).with_context(|| format!("opening {}", path.display()))?)
    };
    let mut raw = BufReader::new(raw);

    let mut inner: Box<dyn BufRead + Send> = if raw.fill_buf()?.starts_with(&GZIP_MAGIC) {
        let decoder = bgzf::io::MultithreadedReader::with_worker_count(threads, raw);
        Box::new(BufReader::new(decoder))
    } else {
        Box::new(raw)
    };

    if inner.fill_buf()?.starts_with(&BCF_MAGIC) {
        let mut reader = bcf::io::Reader::from(inner);
        let header = read_header(|| reader.read_header())?;
        let samples = checked_samples(&header)?;
        let mut record = bcf::Record::default();
        let mut n = 0usize;
        let reader_header = header.clone();
        let records = std::iter::from_fn(move || {
            match reader.read_record(&mut record) {
                Ok(0) => return None,
                Ok(_) => {}
                Err(e) => return Some(Err(e)),
            }
            n += 1;
            let name = record
                .reference_sequence_name(reader_header.string_maps())
                .ok();
            let pos = record.variant_start().and_then(Result::ok).map(usize::from);
            let locus = locus(name, pos, n);
            Some(convert(&reader_header, &record, &locus))
        });
        Ok(VariantSource {
            header,
            samples,
            records: Box::new(records),
        })
    } else {
        let mut reader = vcf::io::Reader::new(inner);
        let header = read_header(|| reader.read_header())?;
        let samples = checked_samples(&header)?;
        let mut record = vcf::Record::default();
        let mut n = 0usize;
        let reader_header = header.clone();
        let records = std::iter::from_fn(move || {
            match reader.read_record(&mut record) {
                Ok(0) => return None,
                Ok(_) => {}
                Err(e) => return Some(Err(e)),
            }
            n += 1;
            let name = Some(record.reference_sequence_name());
            let pos = record.variant_start().and_then(Result::ok).map(usize::from);
            let locus = locus(name, pos, n);
            Some(convert(&reader_header, &record, &locus))
        });
        Ok(VariantSource {
            header,
            samples,
            records: Box::new(records),
        })
    }
}

fn locus(name: Option<&str>, pos: Option<usize>, ordinal: usize) -> String {
    match (name, pos) {
        (Some(name), Some(pos)) => format!("{name}:{pos}"),
        _ => format!("record {ordinal}"),
    }
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

fn convert<R>(header: &Header, record: &R, locus: &str) -> io::Result<RecordBuf>
where
    R: vcf::variant::Record,
{
    let mut buf = RecordBuf::try_from_variant_record(header, record)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("{locus}: {e}")))?;
    decode_record(&mut buf)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("{locus}: {e}")))?;
    Ok(buf)
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

fn decode_all(values: &[Option<String>]) -> io::Result<Vec<Option<String>>> {
    values
        .iter()
        .map(|v| v.as_deref().map(decode).transpose())
        .collect()
}

fn is_missing(array: &SampleArray) -> bool {
    match array {
        SampleArray::Integer(v) => matches!(v.as_slice(), [None]),
        SampleArray::Float(v) => matches!(v.as_slice(), [None]),
        SampleArray::Character(v) => matches!(v.as_slice(), [None]),
        SampleArray::String(v) => matches!(v.as_slice(), [None]),
    }
}

fn decode_record(record: &mut RecordBuf) -> io::Result<()> {
    let info = record.info_mut();
    for i in 0..info.keys().count() {
        if let Some((_, value)) = info.get_index_mut(i) {
            match value {
                Some(InfoValue::String(s)) => *s = decode(s)?,
                Some(InfoValue::Array(InfoArray::String(v))) => *v = decode_all(v)?,
                _ => {}
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
                    Some(SampleValue::Array(SampleArray::String(v))) => Ok(Some(
                        SampleValue::Array(SampleArray::String(decode_all(v)?)),
                    )),
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
    fn opens_every_format(#[values("basic.vcf.gz", "basic.bcf")] name: &str) {
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
