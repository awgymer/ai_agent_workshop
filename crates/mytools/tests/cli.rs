//! Command-line behaviour that needs no bedtools: exit codes for usage errors (SPEC §7)
//! and the determinism of the shared random-BED generator.

mod common;

use common::{Rng, random_bed, run};

const SUBCOMMANDS: [&str; 5] = ["sort", "merge", "intersect", "window", "subtract"];

fn mytools(args: &[&str]) -> common::Output {
    run(common::mytools_path(), args, None)
}

#[test]
fn version_prints_and_exits_0() {
    let out = mytools(&["--version"]);
    assert_eq!(out.code, Some(0));
    assert_eq!(String::from_utf8(out.stdout).unwrap(), format!("mytools {}\n", env!("CARGO_PKG_VERSION")));
}

#[test]
fn no_arguments_is_a_usage_error() {
    let out = mytools(&[]);
    assert_eq!(out.code, Some(2));
    assert!(out.stdout.is_empty());
}

#[test]
fn unknown_subcommand_exits_2_with_empty_stdout() {
    let out = mytools(&["frobnicate"]);
    assert_eq!(out.code, Some(2));
    assert!(out.stdout.is_empty());
}

#[test]
fn unknown_flag_on_every_subcommand_exits_2_with_empty_stdout() {
    for cmd in SUBCOMMANDS {
        let out = mytools(&[cmd, "-zz"]);
        assert_eq!(out.code, Some(2), "{cmd} -zz");
        assert!(out.stdout.is_empty(), "{cmd} -zz");
    }
}

#[test]
fn missing_required_input_exits_2() {
    for cmd in SUBCOMMANDS {
        let out = mytools(&[cmd]);
        assert_eq!(out.code, Some(2), "{cmd} with no arguments");
        assert!(out.stdout.is_empty());
    }
}

#[test]
fn invalid_values_are_rejected_before_input_is_opened() {
    // The input files don't exist: exit 2 rather than 1 proves no input was read.
    let cases: &[&[&str]] = &[
        &["merge", "-i", "nope.bed", "-d", "abc"],
        &["merge", "-i", "nope.bed", "-S", "x"],
        &["merge", "-i", "nope.bed", "-prec", "abc"],
        &["intersect", "-a", "nope.bed", "-b", "nope.bed", "-f", "1.5"],
        &["intersect", "-a", "nope.bed", "-b", "nope.bed", "-f", "-0.5"],
        &["window", "-a", "nope.bed", "-b", "nope.bed", "-w", "abc"],
        &["window", "-a", "nope.bed", "-b", "nope.bed", "-w", "-5"],
        &["subtract", "-a", "nope.bed", "-b", "nope.bed", "-iobuf", "1T"],
    ];
    for args in cases {
        let out = mytools(args);
        assert_eq!(out.code, Some(2), "{args:?}");
        assert!(out.stdout.is_empty(), "{args:?}");
    }
}

#[test]
fn generator_is_deterministic() {
    let a = random_bed(&mut Rng::new(42), 200);
    let b = random_bed(&mut Rng::new(42), 200);
    assert_eq!(a.as_bytes(), b.as_bytes());
    assert_ne!(a, random_bed(&mut Rng::new(43), 200));
}

#[test]
fn generator_makes_valid_edge_case_rich_bed() {
    let bed = random_bed(&mut Rng::new(7), 500);
    let rows: Vec<Vec<&str>> = bed.lines().map(|l| l.split('\t').collect()).collect();
    assert_eq!(rows.len(), 500);
    let mut zero = 0;
    for r in &rows {
        assert_eq!(r.len(), 6);
        assert!(r[0] == "chr1" || r[0] == "chr2");
        let (s, e): (u64, u64) = (r[1].parse().unwrap(), r[2].parse().unwrap());
        assert!(s <= e && e <= common::MAX_COORD, "{r:?}");
        assert!(r[5] == "+" || r[5] == "-");
        zero += usize::from(s == e);
    }
    let coords: Vec<(&str, &str, &str)> = rows.iter().map(|r| (r[0], r[1], r[2])).collect();
    let bookended = coords.iter().filter(|a| coords.iter().any(|b| a.0 == b.0 && a.2 == b.1 && b.1 != b.2)).count();
    let identical = coords.len() - coords.iter().collect::<std::collections::HashSet<_>>().len();
    assert!(zero > 50, "zero-length: {zero}");
    assert!(bookended > 50, "bookended: {bookended}");
    assert!(identical > 20, "identical: {identical}");
}

#[test]
#[ignore = "needs bedtools"]
fn run_both_runs_the_same_arguments() {
    let path = common::temp_file("run_both.bed", &random_bed(&mut Rng::new(1), 20));
    let (_mine, theirs) = common::run_both(&["sort", "-i", path.to_str().unwrap()]);
    assert_eq!(theirs.code, Some(0));
    assert_eq!(theirs.stdout.iter().filter(|&&b| b == b'\n').count(), 20);
}
