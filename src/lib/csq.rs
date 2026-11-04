use anyhow::bail;
use noodles::vcf::Header;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CsqIndex {
    pub names: Vec<String>,
}

impl CsqIndex {
    pub fn from_header(header: &Header) -> Option<CsqIndex> {
        let description = header.infos().get("CSQ")?.description();
        let (_, format) = description.split_once("Format: ")?;
        Some(CsqIndex {
            names: format.trim().split('|').map(str::to_string).collect(),
        })
    }

    pub fn position(&self, name: &str) -> Option<usize> {
        self.names.iter().position(|n| n == name)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CsqEntry {
    pub values: Vec<String>,
}

impl CsqEntry {
    pub fn get<'a>(&'a self, index: &CsqIndex, name: &str) -> Option<&'a str> {
        let value = self.values.get(index.position(name)?)?;
        (!value.is_empty()).then_some(value.as_str())
    }
}

pub fn parse_entries(raw: &[String], index: &CsqIndex) -> anyhow::Result<Vec<CsqEntry>> {
    raw.iter()
        .map(|entry| {
            let values: Vec<String> = entry.split('|').map(str::to_string).collect();
            if values.len() != index.names.len() {
                bail!(
                    "a CSQ entry has {} subfields where the header declares {}",
                    values.len(),
                    index.names.len()
                );
            }
            Ok(CsqEntry { values })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn header(info_lines: &str) -> Header {
        format!("##fileformat=VCFv4.3\n{info_lines}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            .parse()
            .unwrap()
    }

    fn csq_header(description: &str) -> Header {
        header(&format!(
            "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"{description}\">\n"
        ))
    }

    fn index(names: &[&str]) -> CsqIndex {
        CsqIndex {
            names: names.iter().map(ToString::to_string).collect(),
        }
    }

    fn raw(entries: &[&str]) -> Vec<String> {
        entries.iter().map(ToString::to_string).collect()
    }

    #[test]
    fn index_reads_names_after_the_format_marker() {
        let header =
            csq_header("Consequence annotations from VEP. Format: Allele|Consequence|SYMBOL");
        let index = CsqIndex::from_header(&header).unwrap();
        assert_eq!(index.names, ["Allele", "Consequence", "SYMBOL"]);
        assert_eq!(index.position("SYMBOL"), Some(2));
        assert_eq!(index.position("Gene"), None);
    }

    #[test]
    fn index_is_absent_without_a_csq_definition() {
        assert_eq!(CsqIndex::from_header(&header("")), None);
    }

    #[test]
    fn index_is_absent_without_format_text() {
        assert_eq!(
            CsqIndex::from_header(&csq_header("Consequence annotations")),
            None
        );
    }

    #[test]
    fn entries_split_on_pipes_and_empty_is_missing() {
        let index = index(&["Allele", "Consequence", "SYMBOL"]);
        let entries = parse_entries(
            &raw(&["T|missense_variant|", "T|intron_variant|TP53"]),
            &index,
        )
        .unwrap();
        assert_eq!(entries.len(), 2);
        assert_eq!(
            entries[0].get(&index, "Consequence"),
            Some("missense_variant")
        );
        assert_eq!(entries[0].get(&index, "SYMBOL"), None);
        assert_eq!(entries[1].get(&index, "SYMBOL"), Some("TP53"));
        assert_eq!(entries[1].get(&index, "Gene"), None);
    }

    #[test]
    fn entry_with_the_wrong_subfield_count_is_an_error() {
        let index = index(&["Allele", "Consequence", "SYMBOL"]);
        let error = parse_entries(&raw(&["T|a|b|c"]), &index).unwrap_err();
        assert_eq!(
            error.to_string(),
            "a CSQ entry has 4 subfields where the header declares 3"
        );
    }
}
