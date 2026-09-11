//! Sortedness checker for `merge` and the `-sorted` modes (SPEC §4, §7).
//!
//! Input must be sorted by chromosome, then start. bedtools rejects the first record
//! that breaks this, exiting 1; so do we, naming `file:line` and the out-of-order record.
//!
//! What "sorted by chromosome" means depends on the caller (*measured* on `merge`):
//!
//! - [`ChromRule::Grouped`]: each chromosome's records are contiguous, in any chromosome
//!   order. `merge` accepts `chr2` then `chr1`, and `chr10` then `chr2`, but rejects
//!   `chr1`, `chr2`, `chr1`.
//! - [`ChromRule::Ordered`]: chromosomes also follow a [`ChromOrder`], so two files can be
//!   swept together.
//!
//! Memory is one entry per chromosome seen.

use std::cmp::Ordering;
use std::collections::HashSet;

use crate::chrom::ChromOrder;
use crate::error::{Error, Result};
use crate::record::Record;

#[derive(Debug, Clone, Copy)]
pub enum ChromRule<'a> {
    Grouped,
    Ordered(&'a ChromOrder),
}

#[derive(Debug)]
pub struct SortChecker {
    file: String,
    chrom: Option<String>,
    start: u64,
    finished: HashSet<String>,
}

const HINT: &str = "(run mytools sort first)";

impl SortChecker {
    /// `file` is the name used in error messages.
    pub fn new(file: &str) -> SortChecker {
        SortChecker { file: file.to_string(), chrom: None, start: 0, finished: HashSet::new() }
    }

    /// Check the next record; an error names it.
    pub fn check(&mut self, rule: ChromRule, record: &Record) -> Result<()> {
        self.check_position(rule, record.chrom(), record.start, record.line)
    }

    pub fn check_position(&mut self, rule: ChromRule, chrom: &str, start: u64, line: u64) -> Result<()> {
        match &self.chrom {
            Some(prev) if prev == chrom => {
                if start < self.start {
                    return Err(Error::at(
                        &self.file,
                        line,
                        format!("not sorted: {chrom}:{start} comes after {chrom}:{} {HINT}", self.start),
                    ));
                }
            }
            prev => {
                if self.finished.contains(chrom) {
                    let prev = prev.as_deref().unwrap_or_default();
                    return Err(Error::at(
                        &self.file,
                        line,
                        format!("not sorted: {chrom} appears again after {prev} {HINT}"),
                    ));
                }
                if let (ChromRule::Ordered(order), Some(prev)) = (rule, prev) {
                    if order.cmp(prev, chrom) == Ordering::Greater {
                        return Err(Error::at(
                            &self.file,
                            line,
                            format!("not sorted: {chrom} comes after {prev} in chromosome order {HINT}"),
                        ));
                    }
                }
                if let Some(prev) = self.chrom.take() {
                    self.finished.insert(prev);
                }
                self.chrom = Some(chrom.to_string());
            }
        }
        self.start = start;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::Kind;

    fn run(rule: ChromRule, records: &[(&str, u64)]) -> Result<()> {
        let mut c = SortChecker::new("a.bed");
        for (i, &(chrom, start)) in records.iter().enumerate() {
            c.check_position(rule, chrom, start, i as u64 + 1)?;
        }
        Ok(())
    }

    #[test]
    fn sorted_input_passes_including_equal_starts() {
        assert!(run(ChromRule::Grouped, &[("chr1", 0), ("chr1", 500), ("chr1", 500), ("chr2", 0)]).is_ok());
    }

    #[test]
    fn start_going_backwards_names_the_record() {
        let e = run(ChromRule::Grouped, &[("chr1", 300), ("chr1", 0)]).unwrap_err();
        assert_eq!(e.kind(), Kind::Data);
        assert_eq!(e.message(), "a.bed:2: not sorted: chr1:0 comes after chr1:300 (run mytools sort first)");
    }

    #[test]
    fn grouped_accepts_any_chromosome_order() {
        // Measured: bedtools merge accepts chr2 before chr1 and chr10 before chr2.
        assert!(run(ChromRule::Grouped, &[("chr2", 0), ("chr1", 0)]).is_ok());
        assert!(run(ChromRule::Grouped, &[("chr10", 0), ("chr2", 0)]).is_ok());
    }

    #[test]
    fn chromosome_seen_again_is_rejected() {
        let e = run(ChromRule::Grouped, &[("chr1", 0), ("chr2", 0), ("chr1", 20)]).unwrap_err();
        assert_eq!(e.message(), "a.bed:3: not sorted: chr1 appears again after chr2 (run mytools sort first)");
    }

    #[test]
    fn new_chromosome_may_start_below_previous_start() {
        assert!(run(ChromRule::Grouped, &[("chr1", 900), ("chr2", 0)]).is_ok());
    }

    #[test]
    fn ordered_rule_follows_chrom_order() {
        let lex = ChromOrder::Lexicographic;
        assert!(run(ChromRule::Ordered(&lex), &[("chr17", 0), ("chr7", 0)]).is_ok());
        let e = run(ChromRule::Ordered(&lex), &[("chr7", 0), ("chr17", 0)]).unwrap_err();
        assert_eq!(e.message(), "a.bed:2: not sorted: chr17 comes after chr7 in chromosome order (run mytools sort first)");
        let genome = ChromOrder::from_names(["chr7", "chr17"]);
        assert!(run(ChromRule::Ordered(&genome), &[("chr7", 0), ("chr17", 0)]).is_ok());
    }
}
