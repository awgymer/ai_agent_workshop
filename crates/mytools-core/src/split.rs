//! Block splitting for `-split`: BED12 blocks and BAM alignment blocks.
//!
//! Each block is a 0-based, half-open `(start, end)` on the record's chromosome.

use noodles::sam::alignment::record::cigar::op::Kind;

use crate::error::{Error, Result};
use crate::record::Record;

/// The blocks of a record: BAM from its CIGAR, BED12 from columns 10-12, and anything
/// else as one block covering the whole record. `file` names errors.
pub fn blocks(record: &Record, file: &str) -> Result<Vec<(u64, u64)>> {
    if let Some(bam) = &record.bam {
        let mut ops = Vec::new();
        for op in bam.cigar().iter() {
            let op = op.map_err(|e| Error::read(file, &e))?;
            ops.push((op.kind(), op.len() as u64));
        }
        return Ok(cigar_blocks(record.start, &ops, false));
    }
    if record.num_fields() >= 12 {
        return bed12_blocks(record, file);
    }
    Ok(vec![(record.start, record.end)])
}

/// Columns 10-12 of a BED12 record: block count, comma-separated sizes, and starts
/// relative to column 2. A trailing comma is allowed, as UCSC writes one.
pub fn bed12_blocks(record: &Record, file: &str) -> Result<Vec<(u64, u64)>> {
    let bad = |what: String| Error::at(file, record.line, what);
    let field = |i: usize| record.field(i).unwrap_or("");
    let count: usize = field(9)
        .parse()
        .map_err(|_| bad(format!("block count \"{}\" is not an integer", field(9))))?;
    let list = |i: usize, what: &str| -> Result<Vec<u64>> {
        let values: Vec<&str> = field(i).split(',').filter(|v| !v.is_empty()).collect();
        if values.len() != count {
            return Err(bad(format!("expected {count} block {what}, found {}", values.len())));
        }
        values
            .iter()
            .map(|v| v.parse().map_err(|_| bad(format!("block {what} \"{v}\" is not an integer"))))
            .collect()
    };
    let sizes = list(10, "sizes")?;
    let starts = list(11, "starts")?;
    Ok(starts.iter().zip(&sizes).map(|(&s, &len)| (record.start + s, record.start + s + len)).collect())
}

/// Reference blocks covered by a CIGAR starting at 0-based `start`. `M`, `=`, `X` extend
/// a block and `N` ends one. `D` extends the block unless `break_on_deletion`. `I`, `S`,
/// `H` and `P` consume no reference.
pub fn cigar_blocks(start: u64, ops: &[(Kind, u64)], break_on_deletion: bool) -> Vec<(u64, u64)> {
    let mut blocks = Vec::new();
    let mut pos = start;
    let mut block_start = start;
    for &(kind, len) in ops {
        match kind {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => pos += len,
            Kind::Deletion if !break_on_deletion => pos += len,
            Kind::Deletion | Kind::Skip => {
                if pos > block_start {
                    blocks.push((block_start, pos));
                }
                pos += len;
                block_start = pos;
            }
            Kind::Insertion | Kind::SoftClip | Kind::HardClip | Kind::Pad => {}
        }
    }
    if pos > block_start || blocks.is_empty() {
        blocks.push((block_start, pos));
    }
    blocks
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::record::Strand;

    fn bed(text: &str) -> Record {
        let f: Vec<&str> = text.split('\t').collect();
        Record::new(text.to_string(), f[1].parse().unwrap(), f[2].parse().unwrap(), Strand::parse(f[5]), 1)
    }

    #[test]
    fn bed12_blocks_are_absolute() {
        let r = bed("chr1\t100\t1000\ttx1\t0\t+\t150\t900\t0\t3\t100,200,100,\t0,400,800,");
        assert_eq!(blocks(&r, "a.bed12").unwrap(), [(100, 200), (500, 700), (900, 1000)]);
    }

    #[test]
    fn bed12_without_trailing_comma() {
        let r = bed("chr2\t0\t500\ttx3\t0\t+\t0\t500\t0\t2\t100,100\t0,400");
        assert_eq!(blocks(&r, "a.bed12").unwrap(), [(0, 100), (400, 500)]);
    }

    #[test]
    fn bed12_count_mismatch_is_a_data_error() {
        let r = bed("chr1\t0\t500\ttx\t0\t+\t0\t500\t0\t3\t100,100,\t0,400,");
        let e = blocks(&r, "a.bed12").unwrap_err();
        assert_eq!(e.message(), "a.bed12:1: expected 3 block sizes, found 2");
    }

    #[test]
    fn narrower_than_bed12_is_one_block() {
        let r = bed("chr1\t10\t20\tx\t0\t+");
        assert_eq!(blocks(&r, "a.bed").unwrap(), [(10, 20)]);
    }

    #[test]
    fn cigar_skip_splits_and_soft_clip_is_ignored() {
        // 5S10M100N20M
        let ops = [(Kind::SoftClip, 5), (Kind::Match, 10), (Kind::Skip, 100), (Kind::Match, 20)];
        assert_eq!(cigar_blocks(0, &ops, false), [(0, 10), (110, 130)]);
    }

    #[test]
    fn cigar_deletion_extends_unless_asked_to_break() {
        // 10M5D10M, with an insertion that consumes no reference
        let ops = [(Kind::Match, 10), (Kind::Deletion, 5), (Kind::Insertion, 3), (Kind::Match, 10)];
        assert_eq!(cigar_blocks(100, &ops, false), [(100, 125)]);
        assert_eq!(cigar_blocks(100, &ops, true), [(100, 110), (115, 125)]);
    }
}
