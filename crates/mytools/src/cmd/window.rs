//! `mytools window` (#8). Stub: validates the command line, then exits 1.

use mytools_core::args::{self, Flag, Kind, Spec};
use mytools_core::{Error, Result};

const FLAGS: &[Flag] = &[
    Flag::new("a", Kind::Str),
    Flag::new("abam", Kind::Str),
    Flag::new("b", Kind::Str),
    Flag::new("ubam", Kind::Switch),
    Flag::new("bed", Kind::Switch),
    Flag::new("w", Kind::NonNegInt),
    Flag::new("l", Kind::NonNegInt),
    Flag::new("r", Kind::NonNegInt),
    Flag::new("sw", Kind::Switch),
    Flag::new("sm", Kind::Switch),
    Flag::new("Sm", Kind::Switch),
    Flag::new("u", Kind::Switch),
    Flag::new("c", Kind::Switch),
    Flag::new("v", Kind::Switch),
    Flag::new("header", Kind::Switch),
];

const SPEC: Spec = Spec {
    usage: "usage: mytools window [OPTIONS] -a <bed/gff/vcf> | -abam <bam> -b <bed/gff/vcf>",
    flags: FLAGS,
    required: &["a|abam", "b"],
};

pub fn run(argv: &[String]) -> Result<()> {
    let _args = args::parse(&SPEC, argv)?;
    Err(Error::data("not implemented yet"))
}
