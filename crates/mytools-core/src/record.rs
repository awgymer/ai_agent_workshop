//! One interval record, whatever format it was read from.
//!
//! Coordinates are BED's: **0-based, half-open**. GFF and VCF positions are converted on
//! reading (`input.rs`); the original columns are kept untouched in the record's text so
//! they can be carried through to output.

use noodles::bam;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Strand {
    Plus,
    Minus,
    /// `.`, a missing strand column, or anything else.
    Unknown,
}

impl Strand {
    pub fn parse(s: &str) -> Strand {
        match s {
            "+" => Strand::Plus,
            "-" => Strand::Minus,
            _ => Strand::Unknown,
        }
    }

    pub fn as_char(self) -> char {
        match self {
            Strand::Plus => '+',
            Strand::Minus => '-',
            Strand::Unknown => '.',
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Bed,
    Gff,
    Vcf,
    Bam,
}

#[derive(Debug, Clone)]
pub struct Record {
    /// 0-based start.
    pub start: u64,
    /// 0-based, exclusive end. `start == end` is a legal zero-length interval.
    pub end: u64,
    pub strand: Strand,
    /// 1-based line number in the input (record number for BAM).
    pub line: u64,
    /// For BAM input, the alignment itself, so BAM can be written back out.
    pub bam: Option<Box<bam::Record>>,
    /// Tab-separated columns as read. For BAM, the BED6 view of the alignment.
    text: String,
    /// Byte offset just past each column in `text`.
    ends: Vec<u32>,
}

impl Record {
    /// A record whose columns are `text` split on tabs; `text` has no line ending.
    pub fn new(text: String, start: u64, end: u64, strand: Strand, line: u64) -> Record {
        let mut ends: Vec<u32> = text.match_indices('\t').map(|(i, _)| i as u32).collect();
        ends.push(text.len() as u32);
        Record { start, end, strand, line, bam: None, text, ends }
    }

    pub fn chrom(&self) -> &str {
        self.field(0).unwrap_or("")
    }

    /// Column `i`, 0-based.
    pub fn field(&self, i: usize) -> Option<&str> {
        let end = *self.ends.get(i)? as usize;
        let start = if i == 0 { 0 } else { self.ends[i - 1] as usize + 1 };
        Some(&self.text[start..end])
    }

    pub fn num_fields(&self) -> usize {
        self.ends.len()
    }

    pub fn fields(&self) -> impl Iterator<Item = &str> {
        self.text.split('\t')
    }

    /// The columns joined by tabs, exactly as read.
    pub fn text(&self) -> &str {
        &self.text
    }

    pub fn len(&self) -> u64 {
        self.end - self.start
    }

    pub fn is_zero_length(&self) -> bool {
        self.start == self.end
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fields_are_split_on_tabs_and_kept_verbatim() {
        let r = Record::new("chr1\t0\t100\ta01\t10\t+".into(), 0, 100, Strand::Plus, 1);
        assert_eq!(r.chrom(), "chr1");
        assert_eq!(r.field(3), Some("a01"));
        assert_eq!(r.field(5), Some("+"));
        assert_eq!(r.field(6), None);
        assert_eq!(r.num_fields(), 6);
        assert_eq!(r.fields().collect::<Vec<_>>(), ["chr1", "0", "100", "a01", "10", "+"]);
        assert_eq!(r.text(), "chr1\t0\t100\ta01\t10\t+");
    }

    #[test]
    fn empty_columns_are_preserved() {
        let r = Record::new("chr1\t5\t10\t\tx".into(), 5, 10, Strand::Unknown, 1);
        assert_eq!(r.field(3), Some(""));
        assert_eq!(r.field(4), Some("x"));
    }

    #[test]
    fn zero_length_at_position_zero() {
        let r = Record::new("chr2\t0\t0".into(), 0, 0, Strand::Unknown, 1);
        assert!(r.is_zero_length());
        assert_eq!(r.len(), 0);
    }

    #[test]
    fn strand_parsing() {
        assert_eq!(Strand::parse("+"), Strand::Plus);
        assert_eq!(Strand::parse("-"), Strand::Minus);
        assert_eq!(Strand::parse("."), Strand::Unknown);
        assert_eq!(Strand::parse("x"), Strand::Unknown);
    }
}
