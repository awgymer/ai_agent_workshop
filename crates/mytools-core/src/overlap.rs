//! The overlap test, overlap length, `-f`/`-F`/`-r`/`-e` fractions and strand filters
//! (SPEC §3). BED is **0-based, half-open**.
//!
//! Non-zero-length intervals overlap iff `a.start < b.end && b.start < a.end`, strict
//! `<`, so bookended intervals do not overlap. bedtools does not apply that predicate to
//! zero-length intervals as-is: it first widens `[s, s)` to `[s-1, s+1)`. Everything below
//! marked *measured* was observed on bedtools v2.31.1 with `intersect`.

use crate::record::Strand;

/// `[s, s)` becomes `[s-1, s+1)`; anything else is unchanged. Not clamped at 0.
pub fn widen(start: u64, end: u64) -> (i64, i64) {
    if start == end { (start as i64 - 1, end as i64 + 1) } else { (start as i64, end as i64) }
}

/// Whether `a` and `b` overlap, as bedtools decides.
///
/// *Measured* against A = `100-200`: B = `150-150`, `100-100` and `200-200` overlap;
/// `99-99`, `201-201` and bookended `200-300` do not.
pub fn overlaps(a: (u64, u64), b: (u64, u64)) -> bool {
    overlap_bases(a, b) > 0
}

/// Bases shared by the widened intervals; `<= 0` when they don't overlap. This is the
/// numerator of the `-f`/`-F` fractions.
pub fn overlap_bases(a: (u64, u64), b: (u64, u64)) -> i64 {
    let (a_start, a_end) = widen(a.0, a.1);
    let (b_start, b_end) = widen(b.0, b.1);
    a_end.min(b_end) - a_start.max(b_start)
}

/// The overlap `intersect -wo`/`-wao` print for A and B. For ordinary intervals it's the
/// shared bases. For a zero-length B it is the overlap with B widened, minus 2, and A is
/// never widened, so it can be negative. *Measured*, A = `100-200`: B `150-150` -> 0,
/// `100-100` -> -1, `200-200` -> -1; A = `100-100`: B `100-100`, `101-101`, `99-99` -> -2;
/// A `150-150` or `100-100` against B `100-200` -> 0.
pub fn reported_overlap(a: (u64, u64), b: (u64, u64)) -> i64 {
    let (b_start, b_end) = widen(b.0, b.1);
    let shared = (a.1 as i64).min(b_end) - (a.0 as i64).max(b_start);
    if b.0 == b.1 { shared - 2 } else { shared }
}

/// `-f`, `-F`, `-r` and `-e`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Fraction {
    /// `-f`: minimum overlap as a fraction of A.
    pub of_a: f64,
    /// `-F`: minimum overlap as a fraction of B.
    pub of_b: f64,
    /// `-r`: B must also be covered by at least `-f`.
    pub reciprocal: bool,
    /// `-e`: either fraction is enough.
    pub either: bool,
}

impl Default for Fraction {
    /// bedtools' defaults: `1E-9` for both, which is "at least 1 bp" at genomic scales.
    fn default() -> Self {
        Fraction { of_a: 1e-9, of_b: 1e-9, reciprocal: false, either: false }
    }
}

impl Fraction {
    /// Whether `a` and `b` overlap enough. Lengths and overlap use widened coordinates
    /// (*measured*: B = `150-150` in A = `100-200` passes `-F 1.0` and `-f 0.02` but not
    /// `-f 0.021`). bedtools compares in single precision (*measured*: a 1/3 overlap
    /// passes `-f 0.33333335`, which exceeds 1/3 only in double precision), so we do too.
    pub fn accepts(&self, a: (u64, u64), b: (u64, u64)) -> bool {
        let shared = overlap_bases(a, b);
        if shared <= 0 {
            return false;
        }
        let len = |x: (u64, u64)| {
            let (s, e) = widen(x.0, x.1);
            (e - s) as f32
        };
        let a_ok = shared as f32 / len(a) >= self.of_a as f32;
        let b_min = if self.reciprocal { self.of_a } else { self.of_b };
        let b_ok = shared as f32 / len(b) >= b_min as f32;
        if self.either { a_ok || b_ok } else { a_ok && b_ok }
    }
}

/// Strand requirements: `-s`/`-sm` (same), `-S`/`-Sm` (opposite), `merge -S` (one strand).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum StrandFilter {
    #[default]
    Any,
    Same,
    Opposite,
    Only(Strand),
}

impl StrandFilter {
    /// Whether a pair passes. `Same` and `Opposite` never match an unknown strand, not even
    /// `.` against `.` (*measured*, including BED3 against BED3 with `-s`).
    pub fn accepts(self, a: Strand, b: Strand) -> bool {
        let known = a != Strand::Unknown && b != Strand::Unknown;
        match self {
            StrandFilter::Any => true,
            StrandFilter::Same => known && a == b,
            StrandFilter::Opposite => known && a != b,
            StrandFilter::Only(s) => a == s && b == s,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const A: (u64, u64) = (100, 200);

    #[test]
    fn predicate_is_strict_for_ordinary_intervals() {
        assert!(overlaps(A, (150, 250)));
        assert!(overlaps(A, (199, 300)));
        assert!(!overlaps(A, (200, 300)), "bookended after");
        assert!(!overlaps(A, (0, 100)), "bookended before");
        assert!(!overlaps(A, (99, 100)));
        assert!(overlaps(A, (99, 101)));
    }

    #[test]
    fn nested_and_identical_overlap() {
        assert!(overlaps((700, 800), (750, 760)));
        assert!(overlaps((750, 760), (700, 800)));
        assert!(overlaps((700, 800), (700, 800)));
    }

    #[test]
    fn zero_length_b_overlaps_at_both_ends_of_a() {
        assert!(overlaps(A, (150, 150)));
        assert!(overlaps(A, (100, 100)));
        assert!(overlaps(A, (200, 200)));
        assert!(!overlaps(A, (99, 99)));
        assert!(!overlaps(A, (201, 201)));
    }

    #[test]
    fn zero_length_a_is_symmetric() {
        assert!(overlaps((150, 150), A));
        assert!(overlaps((100, 100), A));
        assert!(overlaps((200, 200), A));
        assert!(!overlaps((99, 99), A));
        assert!(!overlaps((201, 201), A));
    }

    #[test]
    fn zero_length_pairs_overlap_within_one_base() {
        assert!(overlaps((100, 100), (100, 100)));
        assert!(overlaps((100, 100), (101, 101)));
        assert!(overlaps((100, 100), (99, 99)));
        assert!(!overlaps((100, 100), (102, 102)));
    }

    #[test]
    fn position_zero() {
        assert!(overlaps((0, 100), (0, 50)));
        assert!(overlaps((0, 0), (0, 5)));
        assert!(!overlaps((0, 0), (1, 5)));
    }

    #[test]
    fn reported_overlap_matches_bedtools_wo() {
        assert_eq!(reported_overlap(A, (150, 250)), 50);
        assert_eq!(reported_overlap(A, (150, 150)), 0);
        assert_eq!(reported_overlap(A, (199, 199)), 0);
        assert_eq!(reported_overlap(A, (101, 101)), 0);
        assert_eq!(reported_overlap(A, (100, 100)), -1);
        assert_eq!(reported_overlap(A, (200, 200)), -1);
        assert_eq!(reported_overlap((150, 150), A), 0);
        assert_eq!(reported_overlap((100, 100), A), 0);
        assert_eq!(reported_overlap((200, 200), A), 0);
        assert_eq!(reported_overlap((100, 100), (100, 100)), -2);
        assert_eq!(reported_overlap((100, 100), (101, 101)), -2);
        assert_eq!(reported_overlap((100, 100), (99, 99)), -2);
    }

    fn frac(of_a: f64, of_b: f64, reciprocal: bool, either: bool) -> Fraction {
        Fraction { of_a, of_b, reciprocal, either }
    }

    #[test]
    fn default_fraction_is_any_overlap() {
        assert!(Fraction::default().accepts(A, (199, 300)));
        assert!(!Fraction::default().accepts(A, (200, 300)));
    }

    #[test]
    fn fraction_of_a_is_inclusive() {
        assert!(frac(0.5, 1e-9, false, false).accepts(A, (150, 200)));
        assert!(!frac(0.51, 1e-9, false, false).accepts(A, (150, 200)));
    }

    #[test]
    fn fraction_uses_single_precision() {
        let f = |of_a| frac(of_a, 1e-9, false, false).accepts((0, 3), (0, 1));
        assert!(f(0.33333335));
        assert!(!f(0.3333334));
    }

    #[test]
    fn fraction_with_zero_length_uses_widened_lengths() {
        assert!(frac(1e-9, 1.0, false, false).accepts(A, (150, 150)));
        assert!(frac(0.02, 1e-9, false, false).accepts(A, (150, 150)));
        assert!(!frac(0.021, 1e-9, false, false).accepts(A, (150, 150)));
        assert!(frac(1e-9, 0.5, false, false).accepts(A, (100, 100)));
        assert!(!frac(1e-9, 0.51, false, false).accepts(A, (100, 100)));
    }

    #[test]
    fn reciprocal_and_either() {
        let b = (190, 400); // 10 bp: 10% of A, ~4.8% of B
        assert!(frac(0.1, 1e-9, false, false).accepts(A, b));
        assert!(!frac(0.1, 1e-9, true, false).accepts(A, b));
        assert!(frac(0.1, 0.9, false, true).accepts(A, b));
        assert!(!frac(0.2, 0.9, false, true).accepts(A, b));
    }

    #[test]
    fn strand_filters_never_match_unknown() {
        use Strand::*;
        assert!(StrandFilter::Any.accepts(Unknown, Minus));
        assert!(StrandFilter::Same.accepts(Plus, Plus));
        assert!(!StrandFilter::Same.accepts(Plus, Minus));
        assert!(!StrandFilter::Same.accepts(Unknown, Unknown));
        assert!(!StrandFilter::Same.accepts(Plus, Unknown));
        assert!(StrandFilter::Opposite.accepts(Plus, Minus));
        assert!(!StrandFilter::Opposite.accepts(Plus, Unknown));
        assert!(!StrandFilter::Opposite.accepts(Unknown, Unknown));
        assert!(StrandFilter::Only(Plus).accepts(Plus, Plus));
        assert!(!StrandFilter::Only(Plus).accepts(Plus, Minus));
    }
}
