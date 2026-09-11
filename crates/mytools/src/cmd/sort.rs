//! `mytools sort` (#5). Stub: validates the command line, then exits 1.

use mytools_core::args::{self, Flag, Kind, Spec};
use mytools_core::{Error, Result};

const FLAGS: &[Flag] = &[
    Flag::new("i", Kind::Str),
    Flag::new("sizeA", Kind::Switch),
    Flag::new("sizeD", Kind::Switch),
    Flag::new("chrThenSizeA", Kind::Switch),
    Flag::new("chrThenSizeD", Kind::Switch),
    Flag::new("chrThenScoreA", Kind::Switch),
    Flag::new("chrThenScoreD", Kind::Switch),
    Flag::new("g", Kind::Str),
    Flag::new("faidx", Kind::Str),
    Flag::new("header", Kind::Switch),
];

const SPEC: Spec = Spec { usage: "usage: mytools sort [OPTIONS] -i <bed/gff/vcf>", flags: FLAGS, required: &["i"] };

pub fn run(argv: &[String]) -> Result<()> {
    let _args = args::parse(&SPEC, argv)?;
    Err(Error::data("not implemented yet"))
}
