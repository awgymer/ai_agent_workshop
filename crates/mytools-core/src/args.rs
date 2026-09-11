//! Hand-written parser for bedtools-style flags (SPEC §4, §9).
//!
//! Flags are a single dash and may be several letters (`-wa`, `-wao`); a value is always
//! the next argument (`-w 500`, `-d -5`). Every value is typed and validated here, so
//! usage errors (exit 2) are raised before any input is opened and stdout stays empty.
//! A flag given twice keeps its last value, as bedtools does.

use std::collections::HashMap;

use crate::error::{Error, Result};

/// What a flag takes.
#[derive(Debug, Clone, Copy)]
pub enum Kind {
    /// No value: present or absent.
    Switch,
    /// Any integer, negative allowed.
    Int,
    /// An integer `>= 0`.
    NonNegInt,
    /// A finite number within a range; each bound may be inclusive or exclusive.
    Float { min: f64, min_inclusive: bool, max: f64, max_inclusive: bool },
    /// Exactly one of the listed strings.
    OneOf(&'static [&'static str]),
    /// A positive byte count with an optional `K`, `M` or `G` suffix (powers of 1024).
    Size,
    /// Any single string (a file name, a column list).
    Str,
    /// One or more strings: every following argument up to the next flag. `-` alone is
    /// a value (stdin), not a flag.
    Many,
}

#[derive(Debug, Clone, Copy)]
pub struct Flag {
    /// Name without the dash: `"wa"` for `-wa`.
    pub name: &'static str,
    pub kind: Kind,
}

impl Flag {
    pub const fn new(name: &'static str, kind: Kind) -> Self {
        Flag { name, kind }
    }
}

/// The command line a subcommand accepts.
#[derive(Debug)]
pub struct Spec {
    /// One usage line, printed after a usage error.
    pub usage: &'static str,
    pub flags: &'static [Flag],
    /// Flags that must be given. `"a|abam"` means at least one of `-a` or `-abam`.
    pub required: &'static [&'static str],
}

#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    Switch,
    Int(i64),
    Float(f64),
    Size(u64),
    Str(String),
    Many(Vec<String>),
}

/// Parsed, validated flags.
#[derive(Debug, Default)]
pub struct Args {
    values: HashMap<&'static str, Value>,
}

impl Args {
    pub fn has(&self, name: &str) -> bool {
        self.values.contains_key(name)
    }

    pub fn int(&self, name: &str) -> Option<i64> {
        match self.values.get(name) {
            Some(Value::Int(v)) => Some(*v),
            _ => None,
        }
    }

    pub fn float(&self, name: &str) -> Option<f64> {
        match self.values.get(name) {
            Some(Value::Float(v)) => Some(*v),
            _ => None,
        }
    }

    pub fn size(&self, name: &str) -> Option<u64> {
        match self.values.get(name) {
            Some(Value::Size(v)) => Some(*v),
            _ => None,
        }
    }

    /// The value of a `Str` or `OneOf` flag, or the first value of a `Many` flag.
    pub fn str(&self, name: &str) -> Option<&str> {
        match self.values.get(name) {
            Some(Value::Str(v)) => Some(v),
            Some(Value::Many(v)) => v.first().map(String::as_str),
            _ => None,
        }
    }

    /// Every value of a `Many` flag; a single-valued flag gives one element.
    pub fn many(&self, name: &str) -> Vec<&str> {
        match self.values.get(name) {
            Some(Value::Many(v)) => v.iter().map(String::as_str).collect(),
            Some(Value::Str(v)) => vec![v.as_str()],
            _ => Vec::new(),
        }
    }
}

/// Parse `argv` (the arguments after the subcommand name) against `spec`.
pub fn parse(spec: &Spec, argv: &[String]) -> Result<Args> {
    let usage = |e: Error| e.with_usage(spec.usage);
    let mut args = Args::default();
    let mut i = 0;
    while i < argv.len() {
        let arg = &argv[i];
        let Some(name) = arg.strip_prefix('-').filter(|n| !n.is_empty()) else {
            return Err(usage(Error::usage(format!("unexpected argument {arg}"))));
        };
        let Some(flag) = spec.flags.iter().find(|f| f.name == name) else {
            return Err(usage(Error::usage(format!("unknown option {arg}"))));
        };
        i += 1;
        let value = match flag.kind {
            Kind::Switch => Value::Switch,
            Kind::Many => {
                let start = i;
                while i < argv.len() && !is_flag(&argv[i]) {
                    i += 1;
                }
                if i == start {
                    return Err(usage(Error::usage(format!("{arg}: missing value"))));
                }
                Value::Many(argv[start..i].to_vec())
            }
            kind => {
                let Some(raw) = argv.get(i) else {
                    return Err(usage(Error::usage(format!("{arg}: missing value"))));
                };
                i += 1;
                parse_value(arg, raw, kind).map_err(usage)?
            }
        };
        args.values.insert(flag.name, value);
    }
    for required in spec.required {
        if !required.split('|').any(|name| args.has(name)) {
            let names: Vec<String> = required.split('|').map(|n| format!("-{n}")).collect();
            return Err(usage(Error::usage(format!("missing required {}", names.join(" or ")))));
        }
    }
    Ok(args)
}

fn is_flag(arg: &str) -> bool {
    arg.len() > 1 && arg.starts_with('-')
}

fn parse_value(flag: &str, raw: &str, kind: Kind) -> Result<Value> {
    let bad = |what: &str| Error::usage(format!("{flag} {raw}: {what}"));
    match kind {
        Kind::Int => raw.parse().map(Value::Int).map_err(|_| bad("must be an integer")),
        Kind::NonNegInt => match raw.parse::<i64>() {
            Ok(v) if v >= 0 => Ok(Value::Int(v)),
            _ => Err(bad("must be a non-negative integer")),
        },
        Kind::Float { min, min_inclusive, max, max_inclusive } => {
            let range = format!(
                "must be a number in {}{min}, {max}{}",
                if min_inclusive { '[' } else { '(' },
                if max_inclusive { ']' } else { ')' },
            );
            match raw.parse::<f64>() {
                Ok(v)
                    if v.is_finite()
                        && (v > min || (min_inclusive && v == min))
                        && (v < max || (max_inclusive && v == max)) =>
                {
                    Ok(Value::Float(v))
                }
                _ => Err(bad(&range)),
            }
        }
        Kind::OneOf(allowed) => {
            if allowed.contains(&raw) {
                Ok(Value::Str(raw.to_string()))
            } else {
                Err(bad(&format!("must be one of {}", allowed.join(" or "))))
            }
        }
        Kind::Size => parse_size(raw)
            .map(Value::Size)
            .ok_or_else(|| bad("must be a positive integer with an optional K, M or G suffix")),
        Kind::Str => Ok(Value::Str(raw.to_string())),
        Kind::Switch | Kind::Many => unreachable!("handled by parse"),
    }
}

/// `4096`, `64K`, `1M`, `2G`. bedtools silently produces no output for `1k`, `0` or
/// `1T`, so those are rejected rather than guessed at.
fn parse_size(raw: &str) -> Option<u64> {
    let (digits, multiplier) = match raw.as_bytes().last()? {
        b'K' => (&raw[..raw.len() - 1], 1 << 10),
        b'M' => (&raw[..raw.len() - 1], 1 << 20),
        b'G' => (&raw[..raw.len() - 1], 1 << 30),
        _ => (raw, 1),
    };
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    let n: u64 = digits.parse().ok()?;
    n.checked_mul(multiplier).filter(|&v| v > 0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::Kind as ErrorKind;

    const FLAGS: &[Flag] = &[
        Flag::new("s", Kind::Switch),
        Flag::new("sorted", Kind::Switch),
        Flag::new("d", Kind::Int),
        Flag::new("w", Kind::NonNegInt),
        Flag::new("f", Kind::Float { min: 0.0, min_inclusive: false, max: 1.0, max_inclusive: true }),
        Flag::new("S", Kind::OneOf(&["+", "-"])),
        Flag::new("iobuf", Kind::Size),
        Flag::new("i", Kind::Str),
        Flag::new("a", Kind::Str),
        Flag::new("abam", Kind::Str),
        Flag::new("b", Kind::Many),
    ];
    const SPEC: Spec = Spec { usage: "usage: test", flags: FLAGS, required: &[] };

    fn argv(s: &str) -> Vec<String> {
        s.split_whitespace().map(String::from).collect()
    }

    fn ok(s: &str) -> Args {
        parse(&SPEC, &argv(s)).unwrap()
    }

    fn usage_err(s: &str) -> String {
        let e = parse(&SPEC, &argv(s)).unwrap_err();
        assert_eq!(e.kind(), ErrorKind::Usage, "{s}");
        assert_eq!(e.usage_line(), Some("usage: test"));
        e.message().to_string()
    }

    #[test]
    fn switch_is_matched_by_exact_name() {
        let a = ok("-sorted");
        assert!(a.has("sorted"));
        assert!(!a.has("s"));
    }

    #[test]
    fn unknown_flag_is_a_usage_error() {
        assert_eq!(usage_err("-zz"), "unknown option -zz");
        assert_eq!(usage_err("stray"), "unexpected argument stray");
    }

    #[test]
    fn int_accepts_negative_values() {
        assert_eq!(ok("-d -5").int("d"), Some(-5));
        assert_eq!(usage_err("-d abc"), "-d abc: must be an integer");
        assert_eq!(usage_err("-d 1.5"), "-d 1.5: must be an integer");
        assert_eq!(usage_err("-d"), "-d: missing value");
    }

    #[test]
    fn non_negative_int_rejects_negative_and_fractions() {
        assert_eq!(ok("-w 0").int("w"), Some(0));
        assert_eq!(usage_err("-w -5"), "-w -5: must be a non-negative integer");
        assert_eq!(usage_err("-w 1.5"), "-w 1.5: must be a non-negative integer");
        assert_eq!(usage_err("-w abc"), "-w abc: must be a non-negative integer");
    }

    #[test]
    fn float_range_respects_open_and_closed_bounds() {
        assert_eq!(ok("-f 1").float("f"), Some(1.0));
        assert_eq!(ok("-f 1E-9").float("f"), Some(1e-9));
        assert_eq!(usage_err("-f 1.5"), "-f 1.5: must be a number in (0, 1]");
        assert_eq!(usage_err("-f 0"), "-f 0: must be a number in (0, 1]");
        assert_eq!(usage_err("-f -0.5"), "-f -0.5: must be a number in (0, 1]");
        assert_eq!(usage_err("-f abc"), "-f abc: must be a number in (0, 1]");
        assert_eq!(usage_err("-f NaN"), "-f NaN: must be a number in (0, 1]");
    }

    #[test]
    fn one_of_accepts_only_listed_values() {
        assert_eq!(ok("-S -").str("S"), Some("-"));
        assert_eq!(ok("-S +").str("S"), Some("+"));
        assert_eq!(usage_err("-S x"), "-S x: must be one of + or -");
    }

    #[test]
    fn size_takes_k_m_g_suffixes() {
        assert_eq!(ok("-iobuf 4096").size("iobuf"), Some(4096));
        assert_eq!(ok("-iobuf 64K").size("iobuf"), Some(64 << 10));
        assert_eq!(ok("-iobuf 1M").size("iobuf"), Some(1 << 20));
        assert_eq!(ok("-iobuf 2G").size("iobuf"), Some(2 << 30));
        for bad in ["1k", "0", "1T", "abc", "1.5M", "-5", "K"] {
            usage_err(&format!("-iobuf {bad}"));
        }
    }

    #[test]
    fn many_collects_until_next_flag_and_keeps_stdin_dash() {
        let a = ok("-b x.bed y.bed - -s");
        assert_eq!(a.many("b"), vec!["x.bed", "y.bed", "-"]);
        assert!(a.has("s"));
        assert_eq!(usage_err("-b -s"), "-b: missing value");
    }

    #[test]
    fn str_value_may_be_a_dash() {
        assert_eq!(ok("-i -").str("i"), Some("-"));
    }

    #[test]
    fn last_occurrence_wins() {
        assert_eq!(ok("-d 1 -d 7").int("d"), Some(7));
    }

    #[test]
    fn required_alternatives() {
        let spec = Spec { usage: "usage: test", flags: FLAGS, required: &["a|abam", "b"] };
        let err = parse(&spec, &argv("-b x")).unwrap_err();
        assert_eq!(err.message(), "missing required -a or -abam");
        assert!(parse(&spec, &argv("-abam x.bam -b y")).is_ok());
        let err = parse(&spec, &argv("-a x")).unwrap_err();
        assert_eq!(err.message(), "missing required -b");
    }
}
