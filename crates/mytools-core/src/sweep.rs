//! One-pass sweep over two sorted inputs, for the `-sorted` modes of `intersect` and
//! `subtract` (SPEC §4, §6).
//!
//! Algorithm: bedtools' chromsweep. `-b` is read in step with `-a`, keeping a cache of B
//! records in file order. For each A record, B records on earlier chromosomes are
//! skipped. The cache is scanned in order: records A starts at or after are dropped,
//! records overlapping A are hits, and the scan stops at the first record starting at or
//! after A's end. Then B records are pulled, tested and cached the same way, stopping at
//! the first that starts at or after A's end. Hits are in B's file order. Starts and ends
//! are widened as in `overlap.rs`.
//!
//! Stopping at the first record past A's end, and dropping by the current A's start, make
//! `-sorted` find fewer hits than the default mode around zero-length intervals. Each of
//! these is *measured*; the default mode reports them all:
//!
//! - A zero-length B behind a longer B with the same start is never reached: A `33-38`
//!   against B `38-47`, `38-38` has no hit.
//! - The same once both are cached: A `37-50`, `47-48` against B `48-50`, `48-48` reports
//!   `48-48` for the first A only.
//! - A zero-length A after an A with the same start misses B ending there: A `49-50`,
//!   `49-49` against B `36-49` has no hit at all.
//!
//! - Time: O(n + m + k) for n A records, m B records and k hits, plus cache scans.
//! - Memory: the B records overlapping the current position, i.e. the overlap depth.
//!
//! Both inputs must be sorted by chromosome in `order` (lexicographic unless `-g`), then
//! by start. B is checked here; the caller checks A with [`SortChecker`] and
//! [`ChromRule::Ordered`] using [`Sweep::order`].
//!
//! Chromosome naming: unless `-nonamecheck`, the first chromosome names seen in A and B
//! must follow the same convention (`chr` prefix or not, zero-padded or not).

use std::cmp::Ordering;
use std::collections::VecDeque;

use crate::chrom::ChromOrder;
use crate::error::{Error, Result};
use crate::input::Input;
use crate::overlap::{overlaps, widen};
use crate::record::Record;
use crate::sorted::{ChromRule, SortChecker};

pub struct Sweep {
    b: Input,
    order: ChromOrder,
    checker: SortChecker,
    next: Option<Record>,
    cache: VecDeque<Record>,
    cache_chrom: String,
    check_names: bool,
}

impl Sweep {
    pub fn new(mut b: Input, order: ChromOrder, nonamecheck: bool) -> Result<Sweep> {
        let mut checker = SortChecker::new(b.name());
        let next = b.next_record()?;
        if let Some(r) = &next {
            checker.check(ChromRule::Ordered(&order), r)?;
        }
        Ok(Sweep {
            b,
            order,
            checker,
            next,
            cache: VecDeque::new(),
            cache_chrom: String::new(),
            check_names: !nonamecheck,
        })
    }

    pub fn order(&self) -> &ChromOrder {
        &self.order
    }

    /// B records overlapping `a`, in B's file order. `a` must not come before the previous
    /// A record.
    pub fn hits(&mut self, a: &Record) -> Result<Vec<&Record>> {
        let chrom = a.chrom();
        if self.check_names {
            if let Some(b) = &self.next {
                check_naming(chrom, b.chrom(), self.b.name())?;
                self.check_names = false;
            }
        }
        if chrom != self.cache_chrom {
            self.cache.clear();
            self.cache_chrom = chrom.to_string();
        }
        while self.next.as_ref().is_some_and(|b| self.order.cmp(b.chrom(), chrom) == Ordering::Less) {
            self.advance()?;
        }
        let a_span = (a.start, a.end);
        let (a_start, a_end) = widen(a.start, a.end);
        let mut hits = Vec::new();
        let mut i = 0;
        while i < self.cache.len() {
            let b_span = (self.cache[i].start, self.cache[i].end);
            let (b_start, b_end) = widen(b_span.0, b_span.1);
            if b_end <= a_start {
                self.cache.remove(i);
                continue;
            }
            if b_start >= a_end {
                break;
            }
            if overlaps(a_span, b_span) {
                hits.push(i);
            }
            i += 1;
        }
        while let Some(b) = &self.next {
            let (b_start, b_end) = widen(b.start, b.end);
            if b.chrom() != chrom || b_start >= a_end {
                break;
            }
            let b = self.next.take().expect("matched Some");
            if b_end > a_start {
                if overlaps(a_span, (b.start, b.end)) {
                    hits.push(self.cache.len());
                }
                self.cache.push_back(b);
            }
            self.advance()?;
        }
        Ok(hits.into_iter().map(|i| &self.cache[i]).collect())
    }

    fn advance(&mut self) -> Result<()> {
        self.next = self.b.next_record()?;
        if let Some(r) = &self.next {
            self.checker.check(ChromRule::Ordered(&self.order), r)?;
        }
        Ok(())
    }
}

/// `chr01` -> (true, true); `1` -> (false, false).
fn convention(chrom: &str) -> (bool, bool) {
    let prefixed = chrom.len() > 3 && chrom[..3].eq_ignore_ascii_case("chr");
    let rest = if prefixed { &chrom[3..] } else { chrom };
    let padded = rest.len() > 1 && rest.starts_with('0') && rest[1..].starts_with(|c: char| c.is_ascii_digit());
    (prefixed, padded)
}

fn check_naming(a: &str, b: &str, b_name: &str) -> Result<()> {
    if convention(a) == convention(b) {
        return Ok(());
    }
    Err(Error::data(format!(
        "{b_name}: chromosome names follow a different convention ({b}) from the other input ({a}); \
         use -nonamecheck to allow"
    )))
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;
    use crate::error::Kind;

    fn input(name: &str, text: &str) -> Input {
        Input::from_reader(name, Box::new(Cursor::new(text.as_bytes().to_vec())), None).unwrap()
    }

    /// For each A record, the 4th column of each B hit.
    fn sweep(a: &str, b: &str, order: ChromOrder) -> Result<Vec<Vec<String>>> {
        let mut a = input("a.bed", a);
        let mut s = Sweep::new(input("b.bed", b), order, false)?;
        let mut out = Vec::new();
        while let Some(rec) = a.next_record()? {
            out.push(s.hits(&rec)?.iter().map(|r| r.field(3).unwrap().to_string()).collect());
        }
        Ok(out)
    }

    #[test]
    fn hits_per_a_record_in_b_order() {
        let a = "chr1\t100\t200\ta1\nchr1\t150\t400\ta2\nchr2\t0\t50\ta3\n";
        let b = "chr1\t0\t1000\tlong\nchr1\t150\t160\tn1\nchr1\t200\t300\tbookend\nchr1\t350\t360\tn2\nchr2\t10\t20\tc2\nchr3\t0\t10\tc3\n";
        let got = sweep(a, b, ChromOrder::Lexicographic).unwrap();
        assert_eq!(got, vec![vec!["long", "n1"], vec!["long", "n1", "bookend", "n2"], vec!["c2"]]);
    }

    #[test]
    fn zero_length_b_at_a_end_is_a_hit_when_it_comes_first() {
        let got = sweep("chr1\t100\t200\ta\n", "chr1\t200\t200\tzero\nchr1\t200\t300\tlong\n", ChromOrder::Lexicographic);
        assert_eq!(got.unwrap(), vec![vec!["zero"]]);
    }

    #[test]
    fn pulling_stops_at_first_b_past_a_end() {
        // bedtools -sorted quirk (measured): the zero-length B behind 38-47 is never
        // reached, though the default mode reports it.
        let got = sweep("chr2\t33\t38\ta\n", "chr2\t38\t47\tlong\nchr2\t38\t38\tzero\n", ChromOrder::Lexicographic);
        assert_eq!(got.unwrap(), vec![Vec::<String>::new()]);
    }

    #[test]
    fn cache_scan_stops_at_first_record_past_a_end() {
        // bedtools -sorted quirk (measured): the second A misses 48-48, which sits in the
        // cache behind 48-50.
        let a = "chr1\t37\t50\ta1\nchr1\t47\t48\ta2\n";
        let got = sweep(a, "chr1\t48\t50\tlong\nchr1\t48\t48\tzero\n", ChromOrder::Lexicographic);
        assert_eq!(got.unwrap(), vec![vec!["long", "zero"], vec![]]);
        // With nothing in front of it, both A records find it (measured).
        let got = sweep(a, "chr1\t48\t48\tzero\n", ChromOrder::Lexicographic);
        assert_eq!(got.unwrap(), vec![vec!["zero"], vec!["zero"]]);
    }

    #[test]
    fn eviction_uses_the_current_a_start() {
        // bedtools -sorted quirk (measured): 49-50 evicts 36-49, so the zero-length A at
        // 49 that follows misses it.
        let got = sweep("chr1\t49\t50\ta1\nchr1\t49\t49\ta2\n", "chr1\t36\t49\tb\n", ChromOrder::Lexicographic);
        assert_eq!(got.unwrap(), vec![Vec::<String>::new(), Vec::new()]);
    }

    #[test]
    fn zero_length_a_sees_b_that_ended_just_before() {
        let a = "chr1\t100\t150\ta1\nchr1\t150\t150\ta2\n";
        let b = "chr1\t120\t150\tb1\n";
        assert_eq!(sweep(a, b, ChromOrder::Lexicographic).unwrap(), vec![vec!["b1"], vec!["b1"]]);
    }

    #[test]
    fn b_chromosomes_missing_from_a_are_skipped() {
        let a = "chr2\t0\t100\ta\n";
        let b = "chr1\t0\t100\tskip\nchr2\t50\t60\thit\n";
        assert_eq!(sweep(a, b, ChromOrder::Lexicographic).unwrap(), vec![vec!["hit"]]);
    }

    #[test]
    fn genome_order_is_followed() {
        let order = ChromOrder::from_names(["chrX", "chr7"]);
        let a = "chrX\t0\t100\ta1\nchr7\t0\t100\ta2\n";
        let b = "chrX\t10\t20\tx\nchr7\t10\t20\tseven\n";
        assert_eq!(sweep(a, b, order).unwrap(), vec![vec!["x"], vec!["seven"]]);
    }

    #[test]
    fn unsorted_b_is_a_data_error() {
        let e = sweep("chr1\t0\t1000\ta\n", "chr1\t500\t600\tb1\nchr1\t100\t200\tb2\n", ChromOrder::Lexicographic)
            .unwrap_err();
        assert_eq!(e.kind(), Kind::Data);
        assert!(e.message().starts_with("b.bed:2: not sorted"), "{}", e.message());
    }

    #[test]
    fn naming_convention_mismatch_needs_nonamecheck() {
        assert!(sweep("chr1\t0\t10\ta\n", "1\t0\t10\tb\n", ChromOrder::Lexicographic).is_err());
        let mut a = input("a.bed", "chr1\t0\t10\ta\n");
        let mut s = Sweep::new(input("b.bed", "1\t0\t10\tb\n"), ChromOrder::Lexicographic, true).unwrap();
        let rec = a.next_record().unwrap().unwrap();
        assert!(s.hits(&rec).unwrap().is_empty());
        assert_eq!(convention("chr01"), (true, true));
        assert_eq!(convention("chr10"), (true, false));
        assert_eq!(convention("X"), (false, false));
    }
}
