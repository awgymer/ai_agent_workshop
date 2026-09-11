//! Chromosome ordering: lexicographic, or the order listed in a genome file (`-g`) or a
//! FASTA index (`-faidx`).
//!
//! Lexicographic means byte order, as bedtools compares `std::string`s: `chr17` sorts
//! before `chr7`, and `chrX` after both.

use std::cmp::Ordering;
use std::collections::HashMap;
use std::fs;

use crate::error::{Error, Result};

#[derive(Debug, Clone, Default)]
pub enum ChromOrder {
    #[default]
    Lexicographic,
    /// Chromosomes in the order a file lists them.
    Listed(HashMap<String, usize>),
}

impl ChromOrder {
    /// Read a genome file (`name<TAB>size`) or `.fai`: the first column of each non-empty
    /// line gives the order. A name listed twice keeps its first position.
    pub fn from_file(path: &str) -> Result<ChromOrder> {
        let text = fs::read_to_string(path).map_err(|e| Error::open(path, &e))?;
        Ok(ChromOrder::from_names(
            text.lines().filter_map(|l| l.split('\t').next()).filter(|n| !n.trim().is_empty()),
        ))
    }

    pub fn from_names<'a>(names: impl IntoIterator<Item = &'a str>) -> ChromOrder {
        let mut ranks = HashMap::new();
        for name in names {
            let next = ranks.len();
            ranks.entry(name.to_string()).or_insert(next);
        }
        ChromOrder::Listed(ranks)
    }

    /// Position in the listed order, or `None` for lexicographic order and for names the
    /// file doesn't list.
    pub fn rank(&self, chrom: &str) -> Option<usize> {
        match self {
            ChromOrder::Lexicographic => None,
            ChromOrder::Listed(ranks) => ranks.get(chrom).copied(),
        }
    }

    /// Whether `chrom` has a place in this order. Always true for lexicographic order.
    pub fn contains(&self, chrom: &str) -> bool {
        match self {
            ChromOrder::Lexicographic => true,
            ChromOrder::Listed(ranks) => ranks.contains_key(chrom),
        }
    }

    /// Compare two chromosome names. In a listed order, names the file doesn't list sort
    /// after every listed name, lexicographically among themselves.
    pub fn cmp(&self, a: &str, b: &str) -> Ordering {
        match self {
            ChromOrder::Lexicographic => a.as_bytes().cmp(b.as_bytes()),
            ChromOrder::Listed(ranks) => match (ranks.get(a), ranks.get(b)) {
                (Some(x), Some(y)) => x.cmp(y),
                (Some(_), None) => Ordering::Less,
                (None, Some(_)) => Ordering::Greater,
                (None, None) => a.as_bytes().cmp(b.as_bytes()),
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn genome_and_fai_fixtures_give_the_same_order() {
        let dir = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/fixtures/");
        for file in ["genome.txt", "genome.fai"] {
            let o = ChromOrder::from_file(&format!("{dir}{file}")).unwrap();
            let ranks: Vec<_> = ["chr1", "chr2", "chr3", "chrX", "chr7", "chr17"].iter().map(|c| o.rank(c)).collect();
            assert_eq!(ranks, [Some(0), Some(1), Some(2), Some(3), Some(4), Some(5)], "{file}");
        }
        let e = ChromOrder::from_file("nope.genome").unwrap_err();
        assert_eq!(e.message(), "cannot open nope.genome: No such file or directory");
    }

    #[test]
    fn lexicographic_puts_chr17_before_chr7() {
        let o = ChromOrder::Lexicographic;
        assert_eq!(o.cmp("chr17", "chr7"), Ordering::Less);
        assert_eq!(o.cmp("chr7", "chrX"), Ordering::Less);
        assert_eq!(o.cmp("chr1", "chr1"), Ordering::Equal);
        // Byte order: upper case before lower case.
        assert_eq!(o.cmp("chrX", "chr_alt"), Ordering::Less);
    }

    #[test]
    fn listed_order_follows_the_file() {
        let o = ChromOrder::from_names(["chr1", "chr2", "chr3", "chrX", "chr7", "chr17"]);
        assert_eq!(o.cmp("chrX", "chr7"), Ordering::Less);
        assert_eq!(o.cmp("chr7", "chr17"), Ordering::Less);
        assert_eq!(o.rank("chr3"), Some(2));
        assert!(o.contains("chr17"));
        assert!(!o.contains("chrM"));
    }

    #[test]
    fn unlisted_names_sort_last() {
        let o = ChromOrder::from_names(["chr2", "chr1"]);
        assert_eq!(o.cmp("chrM", "chr1"), Ordering::Greater);
        assert_eq!(o.cmp("chrA", "chrM"), Ordering::Less);
        assert_eq!(o.rank("chrM"), None);
    }

    #[test]
    fn duplicate_name_keeps_first_position() {
        let o = ChromOrder::from_names(["chr2", "chr1", "chr2"]);
        assert_eq!(o.rank("chr2"), Some(0));
        assert_eq!(o.rank("chr1"), Some(1));
    }
}
