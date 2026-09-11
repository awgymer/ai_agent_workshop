//! Everything the `mytools` subcommands share: record I/O, argument parsing, errors and
//! exit codes, chromosome ordering, the overlap test, the interval index, the sorted
//! sweep and interval merging. `SPEC.md` is the contract; BED is 0-based, half-open.

pub mod args;
pub mod chrom;
pub mod error;
pub mod index;
pub mod input;
pub mod merge;
pub mod output;
pub mod overlap;
pub mod record;
pub mod sorted;
pub mod split;
pub mod sweep;

pub use error::{Error, Result};
