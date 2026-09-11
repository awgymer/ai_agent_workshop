//! `mytools merge` (#6). Stub: validates the command line, then exits 1.

use mytools_core::args::{self, Flag, Kind, Spec};
use mytools_core::{Error, Result};

const FLAGS: &[Flag] = &[
    Flag::new("i", Kind::Str),
    Flag::new("s", Kind::Switch),
    Flag::new("S", Kind::OneOf(&["+", "-"])),
    Flag::new("d", Kind::Int),
    Flag::new("c", Kind::Str),
    Flag::new("o", Kind::Str),
    Flag::new("delim", Kind::Str),
    Flag::new("prec", Kind::NonNegInt),
    Flag::new("bed", Kind::Switch),
    Flag::new("header", Kind::Switch),
    Flag::new("nobuf", Kind::Switch),
    Flag::new("iobuf", Kind::Size),
];

const SPEC: Spec = Spec { usage: "usage: mytools merge [OPTIONS] -i <bed/gff/vcf/bam>", flags: FLAGS, required: &["i"] };

pub fn run(argv: &[String]) -> Result<()> {
    let _args = args::parse(&SPEC, argv)?;
    Err(Error::data("not implemented yet"))
}
