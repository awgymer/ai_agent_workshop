//! Per-chromosome binned interval index, for the default (unsorted) modes of `intersect`,
//! `window` and `subtract` (SPEC §6).
//!
//! Algorithm: bedtools' hierarchical binning. Level 0 bins are 16 kb (`>> 14`) and each
//! level up is 8 times larger (`>> 3`), seven levels in all, the last covering 4 Gb. A
//! record lives in the smallest bin that holds it whole. A query visits, level by level,
//! every bin its interval touches.
//!
//! - Build: O(m) time and O(m) memory for m records.
//! - Query: time proportional to the records in the bins it touches.
//!
//! Hits come back in bedtools' order (*measured* with `intersect -wb`): finest level
//! first, lower bins first within a level, and file order within a bin. Zero-length
//! records are binned and matched with the widening described in `overlap.rs`.
//!
//! A zero-length record at position 0 widens to `[-1, 1)`, which bedtools' index can't
//! bin: `intersect` and `subtract` exit 1 ("illegal bin number -1") when `-b` has one,
//! while `window` and the `-sorted` modes accept it (*measured*). [`IntervalIndex::from_input`]
//! rejects it for the former; `window` should build with [`IntervalIndex::insert`].

use std::collections::HashMap;

use crate::error::{Error, Result};
use crate::input::Input;
use crate::overlap::{overlaps, widen};
use crate::record::Record;

const FIRST_SHIFT: u32 = 14;
const NEXT_SHIFT: u32 = 3;
const LEVELS: u32 = 7;

#[derive(Debug, Default)]
pub struct IntervalIndex {
    records: Vec<Record>,
    /// chromosome -> (level, bin) key -> record indices in insertion order.
    bins: HashMap<String, HashMap<u64, Vec<u32>>>,
}

fn key(level: u32, bin: u64) -> u64 {
    ((level as u64) << 56) | bin
}

/// The widened interval as bin coordinates: `[first, last]` positions, clamped at 0.
fn span(start: u64, end: u64) -> (u64, u64) {
    let (s, e) = widen(start, end);
    let first = s.max(0) as u64;
    let last = (e - 1).max(0) as u64;
    (first, last.max(first))
}

impl IntervalIndex {
    pub fn new() -> IntervalIndex {
        IntervalIndex::default()
    }

    /// Read every record of `input` into an index, rejecting a zero-length record at
    /// position 0 as `intersect` and `subtract` do (see the module notes).
    pub fn from_input(input: &mut Input) -> Result<IntervalIndex> {
        let mut index = IntervalIndex::new();
        while let Some(record) = input.next_record()? {
            if record.start == 0 && record.end == 0 {
                return Err(Error::at(
                    input.name(),
                    record.line,
                    format!("zero-length interval at {}:0 can't be indexed (bedtools rejects it too)", record.chrom()),
                ));
            }
            index.insert(record);
        }
        Ok(index)
    }

    /// Add one record. Never fails: a zero-length record at 0 is binned at 0.
    pub fn insert(&mut self, record: Record) {
        let (first, last) = span(record.start, record.end);
        let mut level = LEVELS - 1;
        let mut shift = FIRST_SHIFT;
        for l in 0..LEVELS {
            if first >> shift == last >> shift {
                level = l;
                break;
            }
            shift += NEXT_SHIFT;
        }
        let shift = FIRST_SHIFT + NEXT_SHIFT * level;
        let id = self.records.len() as u32;
        self.bins
            .entry(record.chrom().to_string())
            .or_default()
            .entry(key(level, first >> shift))
            .or_default()
            .push(id);
        self.records.push(record);
    }

    pub fn len(&self) -> usize {
        self.records.len()
    }

    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// Every record, in insertion order.
    pub fn records(&self) -> &[Record] {
        &self.records
    }

    /// Records on `chrom` that overlap `[start, end)`, in bedtools' order.
    pub fn query(&self, chrom: &str, start: u64, end: u64) -> Vec<&Record> {
        let mut hits = Vec::new();
        self.for_each_hit(chrom, start, end, |r| hits.push(r));
        hits
    }

    /// As [`IntervalIndex::query`], calling `f` for each hit instead of collecting them.
    pub fn for_each_hit<'a>(&'a self, chrom: &str, start: u64, end: u64, mut f: impl FnMut(&'a Record)) {
        let Some(bins) = self.bins.get(chrom) else { return };
        let (first, last) = span(start, end);
        let mut shift = FIRST_SHIFT;
        for level in 0..LEVELS {
            for bin in (first >> shift)..=(last >> shift) {
                for &id in bins.get(&key(level, bin)).into_iter().flatten() {
                    let r = &self.records[id as usize];
                    if overlaps((start, end), (r.start, r.end)) {
                        f(r);
                    }
                }
            }
            shift += NEXT_SHIFT;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::record::Strand;

    fn index(records: &[(&str, u64, u64, &str)]) -> IntervalIndex {
        let mut idx = IntervalIndex::new();
        for (i, &(chrom, s, e, name)) in records.iter().enumerate() {
            idx.insert(Record::new(format!("{chrom}\t{s}\t{e}\t{name}"), s, e, Strand::Unknown, i as u64 + 1));
        }
        idx
    }

    fn names(idx: &IntervalIndex, chrom: &str, s: u64, e: u64) -> Vec<String> {
        idx.query(chrom, s, e).iter().map(|r| r.field(3).unwrap().to_string()).collect()
    }

    const BIG: (u64, u64) = (0, 100_000_000);

    #[test]
    fn hit_order_matches_bedtools() {
        // Each case measured with `bedtools intersect -wb -a big.bed -b <these, in order>`.
        let cases: &[(&[(&str, u64, u64, &str)], &[&str])] = &[
            (&[("chr1", 10, 16385, "cross"), ("chr1", 10, 16384, "inside")], &["inside", "cross"]),
            (&[("chr1", 10, 131073, "cross"), ("chr1", 10, 131072, "inside")], &["inside", "cross"]),
            (&[("chr1", 500, 600, "x"), ("chr1", 10, 20, "y")], &["x", "y"]),
            (&[("chr1", 20000, 20010, "hi"), ("chr1", 10, 20, "lo")], &["lo", "hi"]),
            (
                &[("chr1", 10, 2_000_000, "L2"), ("chr1", 10, 200_000, "L1"), ("chr1", 10, 20, "L0")],
                &["L0", "L1", "L2"],
            ),
            (
                &[("chr1", 500, 600, "b1"), ("chr1", 10, 90000, "b2"), ("chr1", 100, 200, "b3"), ("chr1", 50, 60, "b4")],
                &["b1", "b3", "b4", "b2"],
            ),
        ];
        for (records, expected) in cases {
            assert_eq!(names(&index(records), "chr1", BIG.0, BIG.1), *expected);
        }
    }

    #[test]
    fn query_filters_with_the_overlap_predicate() {
        let idx = index(&[
            ("chr1", 200, 300, "bookended"),
            ("chr1", 150, 150, "zero"),
            ("chr1", 200, 200, "zero_at_end"),
            ("chr1", 0, 100, "before"),
            ("chr1", 120, 130, "nested"),
            ("chr2", 150, 160, "other_chrom"),
        ]);
        assert_eq!(names(&idx, "chr1", 100, 200), ["zero", "zero_at_end", "nested"]);
        assert!(names(&idx, "chr3", 0, 1000).is_empty());
    }

    #[test]
    fn position_zero_and_zero_length_query() {
        let idx = index(&[("chr1", 0, 5, "p"), ("chr1", 1, 5, "q")]);
        assert_eq!(names(&idx, "chr1", 0, 0), ["p"]);
    }

    #[test]
    fn from_input_rejects_zero_length_at_position_zero() {
        let text = "chr1\t0\t10\nchr2\t0\t0\n";
        let mut input = Input::from_reader("b.bed", Box::new(std::io::Cursor::new(text.as_bytes().to_vec())), None).unwrap();
        let e = IntervalIndex::from_input(&mut input).unwrap_err();
        assert_eq!(e.exit_code(), 1);
        assert!(e.message().starts_with("b.bed:2: zero-length interval at chr2:0"), "{}", e.message());
        let mut ok = Input::from_reader("b.bed", Box::new(std::io::Cursor::new(b"chr1\t1\t1\n".to_vec())), None).unwrap();
        assert_eq!(IntervalIndex::from_input(&mut ok).unwrap().len(), 1);
    }

    #[test]
    fn query_spanning_level0_bins_finds_records_on_both_sides() {
        let idx = index(&[("chr1", 16000, 16100, "left"), ("chr1", 16400, 16500, "right")]);
        assert_eq!(names(&idx, "chr1", 16050, 16450), ["left", "right"]);
    }
}
