//! Helpers shared by the integration tests: a seeded random BED generator and
//! [`run_both`], which runs the built `mytools` and real bedtools on the same arguments.
//!
//! Use from a test file with `mod common;`.

#![allow(dead_code)]

use std::fmt::Write as _;
use std::io::Write as _;
use std::path::PathBuf;
use std::process::{Command, Stdio};

/// SplitMix64. Tiny and fully specified, so a seed gives the same inputs everywhere.
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Rng {
        Rng(seed)
    }

    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    /// Uniform in `0..n`.
    pub fn below(&mut self, n: u64) -> u64 {
        self.next_u64() % n
    }

    pub fn pick<'a, T>(&mut self, items: &'a [T]) -> &'a T {
        &items[self.below(items.len() as u64) as usize]
    }
}

/// Largest coordinate the generator uses. Packing everything into 0-50 makes overlaps,
/// bookends and duplicates common.
pub const MAX_COORD: u64 = 50;

/// `n` random BED6 records on `chr1`/`chr2`, coordinates in `0..=50`, both strands, with
/// plenty of zero-length, bookended, nested and identical intervals. Unsorted: sort with
/// `bedtools sort` for commands that need it.
pub fn random_bed(rng: &mut Rng, n: usize) -> String {
    let mut records: Vec<(&str, u64, u64)> = Vec::with_capacity(n);
    let mut out = String::new();
    for i in 0..n {
        let previous = if records.is_empty() { None } else { Some(*rng.pick(&records)) };
        let (chrom, start, end) = match (rng.below(10), previous) {
            (0 | 1, _) => {
                let s = rng.below(MAX_COORD + 1);
                (*rng.pick(&["chr1", "chr2"]), s, s)
            }
            (2, Some((c, _, e))) => (c, e, (e + 1 + rng.below(10)).min(MAX_COORD)),
            (3, Some((c, s, _))) => (c, s.saturating_sub(1 + rng.below(10)), s),
            (4 | 5, Some((c, s, e))) => {
                let ns = s + rng.below(e - s + 1);
                (c, ns, ns + rng.below(e - ns + 1))
            }
            (6, Some(p)) => p,
            _ => {
                let s = rng.below(MAX_COORD);
                (*rng.pick(&["chr1", "chr2"]), s, s + rng.below(MAX_COORD - s + 1))
            }
        };
        records.push((chrom, start, end));
        let score = rng.below(100);
        let strand = rng.pick(&["+", "-"]);
        writeln!(out, "{chrom}\t{start}\t{end}\tr{i}\t{score}\t{strand}").unwrap();
    }
    out
}

/// What a finished process produced.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Output {
    pub stdout: Vec<u8>,
    /// `None` if killed by a signal.
    pub code: Option<i32>,
}

pub fn mytools_path() -> &'static str {
    env!("CARGO_BIN_EXE_mytools")
}

/// Run `program args...`, optionally feeding `stdin`.
pub fn run(program: &str, args: &[&str], stdin: Option<&[u8]>) -> Output {
    let mut child = Command::new(program)
        .args(args)
        .stdin(if stdin.is_some() { Stdio::piped() } else { Stdio::null() })
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap_or_else(|e| panic!("cannot run {program}: {e}"));
    if let Some(bytes) = stdin {
        child.stdin.take().unwrap().write_all(bytes).unwrap();
    }
    let out = child.wait_with_output().unwrap();
    Output { stdout: out.stdout, code: out.status.code() }
}

/// Run the built `mytools` and `bedtools` with the same arguments: `(mytools, bedtools)`.
pub fn run_both(args: &[&str]) -> (Output, Output) {
    (run(mytools_path(), args, None), run("bedtools", args, None))
}

/// Write `contents` to a file under Cargo's per-target temp directory and return its path.
pub fn temp_file(name: &str, contents: &str) -> PathBuf {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join(format!("mytools-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join(name);
    std::fs::write(&path, contents).unwrap();
    path
}
