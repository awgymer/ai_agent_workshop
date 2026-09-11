//! Merging sorted intervals within distance `d` (SPEC §3).
//!
//! Algorithm: one pass. The current group is extended while the next interval starts no
//! more than `d` bases after the group's end: `next.start - group.end <= d`. So `d = 0`
//! merges overlapping *and* bookended intervals, and a negative `d` requires that many
//! bases of overlap. O(n) time, O(1) memory.
//!
//! Zero-length intervals follow bedtools v2.31.1, which widens `[s, s)` to `[s-1, s+1)`
//! before comparing (*measured*):
//!
//! - `0-0` + `0-5` gives `-1 5`: widening is not clamped at position 0.
//! - `5-5` + `8-10` stays apart at `-d 1` but gives `4 10` at `-d 2`.
//! - A group of one record prints its original coordinates: `chr2 0 0` stays `chr2 0 0`.
//!
//! The caller checks sortedness and decides what "same group" means for strand.

use std::cmp::max;

use crate::overlap::widen;

/// A finished group. `start` can be `-1` (see the module notes).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Merged {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    /// Number of input intervals in the group.
    pub count: usize,
}

#[derive(Debug)]
struct Group {
    chrom: String,
    start: i64,
    end: i64,
    original: (u64, u64),
    count: usize,
}

impl Group {
    fn finish(self) -> Merged {
        let (start, end) = if self.count == 1 {
            (self.original.0 as i64, self.original.1 as i64)
        } else {
            (self.start, self.end)
        };
        Merged { chrom: self.chrom, start, end, count: self.count }
    }
}

#[derive(Debug)]
pub struct Merger {
    d: i64,
    group: Option<Group>,
}

impl Merger {
    pub fn new(d: i64) -> Merger {
        Merger { d, group: None }
    }

    /// Add the next interval. Returns the group it closed, if it didn't join the current
    /// one; the interval then starts a new group.
    pub fn push(&mut self, chrom: &str, start: u64, end: u64) -> Option<Merged> {
        let (ws, we) = widen(start, end);
        if let Some(g) = &mut self.group {
            if g.chrom == chrom && ws - g.end <= self.d {
                // The group keeps its first record's (widened) start; a zero-length record
                // joining later doesn't pull it left. The end is the furthest seen.
                g.end = max(g.end, we);
                g.count += 1;
                return None;
            }
        }
        let new = Group { chrom: chrom.to_string(), start: ws, end: we, original: (start, end), count: 1 };
        self.group.replace(new).map(Group::finish)
    }

    /// Close the last group.
    pub fn finish(&mut self) -> Option<Merged> {
        self.group.take().map(Group::finish)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn merge(d: i64, intervals: &[(&str, u64, u64)]) -> Vec<(String, i64, i64, usize)> {
        let mut m = Merger::new(d);
        let mut out = Vec::new();
        for &(c, s, e) in intervals {
            out.extend(m.push(c, s, e));
        }
        out.extend(m.finish());
        out.into_iter().map(|g| (g.chrom, g.start, g.end, g.count)).collect()
    }

    fn coords(d: i64, intervals: &[(u64, u64)]) -> Vec<(i64, i64)> {
        let with_chrom: Vec<_> = intervals.iter().map(|&(s, e)| ("chr1", s, e)).collect();
        merge(d, &with_chrom).into_iter().map(|(_, s, e, _)| (s, e)).collect()
    }

    #[test]
    fn bookended_merge_at_d0() {
        assert_eq!(coords(0, &[(100, 200), (200, 300)]), [(100, 300)]);
    }

    #[test]
    fn gap_of_one_not_merged_at_d0() {
        assert_eq!(coords(0, &[(100, 200), (201, 300)]), [(100, 200), (201, 300)]);
        assert_eq!(coords(1, &[(100, 200), (201, 300)]), [(100, 300)]);
    }

    #[test]
    fn negative_d_requires_overlap() {
        // Measured on bedtools: 0-10, 8-20, 20-30, 25-40.
        let input = [(0, 10), (8, 20), (20, 30), (25, 40)];
        assert_eq!(coords(0, &input), [(0, 40)]);
        assert_eq!(coords(-1, &input), [(0, 20), (20, 40)]);
        assert_eq!(coords(-2, &input), [(0, 20), (20, 40)]);
        assert_eq!(coords(-3, &input), [(0, 10), (8, 20), (20, 40)]);
        assert_eq!(coords(-5, &input), [(0, 10), (8, 20), (20, 40)]);
    }

    #[test]
    fn nested_and_identical() {
        assert_eq!(coords(0, &[(300, 400), (320, 350)]), [(300, 400)]);
        assert_eq!(merge(0, &[("chr1", 700, 800), ("chr1", 700, 800)]), [("chr1".into(), 700, 800, 2)]);
    }

    #[test]
    fn position_zero() {
        assert_eq!(coords(0, &[(0, 100), (100, 200)]), [(0, 200)]);
    }

    #[test]
    fn chromosome_change_closes_group() {
        let out = merge(0, &[("chr1", 0, 100), ("chr2", 50, 150)]);
        assert_eq!(out, [("chr1".into(), 0, 100, 1), ("chr2".into(), 50, 150, 1)]);
    }

    // Zero-length intervals: bedtools widens them by one base each side. Every expected
    // value below was measured on bedtools v2.31.1; none of it is what you'd guess.
    #[test]
    fn zero_length_widens_left_into_neighbour() {
        // data/a.bed: a07 chr1 500 500 + a08 chr1 500 600 -> chr1 499 600.
        assert_eq!(coords(0, &[(500, 500), (500, 600)]), [(499, 600)]);
    }

    #[test]
    fn lone_zero_length_keeps_original_coordinates() {
        assert_eq!(coords(0, &[(0, 0)]), [(0, 0)]);
        assert_eq!(coords(0, &[(1, 1)]), [(1, 1)]);
        assert_eq!(coords(-1, &[(5, 5), (6, 10)]), [(5, 5), (6, 10)]);
    }

    #[test]
    fn zero_length_widening_is_not_clamped_at_zero() {
        assert_eq!(coords(0, &[(0, 0), (0, 5)]), [(-1, 5)]);
        assert_eq!(coords(0, &[(0, 0), (1, 5)]), [(-1, 5)]);
        assert_eq!(coords(0, &[(0, 0), (0, 0)]), [(-1, 1)]);
    }

    #[test]
    fn zero_length_widens_right_past_end() {
        assert_eq!(coords(0, &[(0, 10), (10, 10)]), [(0, 11)]);
        assert_eq!(coords(0, &[(0, 10), (11, 11)]), [(0, 12)]);
        assert_eq!(coords(0, &[(0, 10), (5, 5)]), [(0, 10)]);
        assert_eq!(coords(0, &[(3, 3), (3, 3)]), [(2, 4)]);
    }

    #[test]
    fn zero_length_joining_later_does_not_lower_group_start() {
        // Found by a random comparison: bedtools gives 49 50, not 48 50.
        assert_eq!(coords(-2, &[(49, 50), (49, 49)]), [(49, 50)]);
        assert_eq!(coords(-2, &[(49, 49), (49, 50)]), [(49, 49), (49, 50)]);
        assert_eq!(coords(-2, &[(48, 50), (49, 49)]), [(48, 50)]);
    }

    #[test]
    fn zero_length_widening_counts_towards_distance() {
        assert_eq!(coords(1, &[(5, 5), (8, 10)]), [(5, 5), (8, 10)]);
        assert_eq!(coords(2, &[(5, 5), (8, 10)]), [(4, 10)]);
        assert_eq!(coords(-1, &[(20, 20), (21, 21)]), [(19, 22)]);
    }
}
