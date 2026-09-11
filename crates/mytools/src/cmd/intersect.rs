//! `mytools intersect` (#7). Stub: validates the command line, then exits 1.

use mytools_core::args::{self, Flag, Kind, Spec};
use mytools_core::{Error, Result};

const FRACTION: Kind = Kind::Float { min: 0.0, min_inclusive: false, max: 1.0, max_inclusive: true };

const FLAGS: &[Flag] = &[
    Flag::new("a", Kind::Str),
    Flag::new("b", Kind::Many),
    Flag::new("wa", Kind::Switch),
    Flag::new("wb", Kind::Switch),
    Flag::new("loj", Kind::Switch),
    Flag::new("wo", Kind::Switch),
    Flag::new("wao", Kind::Switch),
    Flag::new("u", Kind::Switch),
    Flag::new("c", Kind::Switch),
    Flag::new("C", Kind::Switch),
    Flag::new("v", Kind::Switch),
    Flag::new("ubam", Kind::Switch),
    Flag::new("s", Kind::Switch),
    Flag::new("S", Kind::Switch),
    Flag::new("f", FRACTION),
    Flag::new("F", FRACTION),
    Flag::new("r", Kind::Switch),
    Flag::new("e", Kind::Switch),
    Flag::new("split", Kind::Switch),
    Flag::new("g", Kind::Str),
    Flag::new("nonamecheck", Kind::Switch),
    Flag::new("sorted", Kind::Switch),
    Flag::new("names", Kind::Many),
    Flag::new("filenames", Kind::Switch),
    Flag::new("sortout", Kind::Switch),
    Flag::new("bed", Kind::Switch),
    Flag::new("header", Kind::Switch),
    Flag::new("nobuf", Kind::Switch),
    Flag::new("iobuf", Kind::Size),
];

const SPEC: Spec = Spec {
    usage: "usage: mytools intersect [OPTIONS] -a <bed/gff/vcf/bam> -b <bed/gff/vcf/bam>...",
    flags: FLAGS,
    required: &["a", "b"],
};

pub fn run(argv: &[String]) -> Result<()> {
    let _args = args::parse(&SPEC, argv)?;
    Err(Error::data("not implemented yet"))
}
