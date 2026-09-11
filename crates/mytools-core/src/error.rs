//! The error type every subcommand returns, carrying its exit code (SPEC §7).
//!
//! - exit `1`: bad input data. The command line was valid, but an input can't be opened
//!   or parsed. The message names `file:line` where there is one.
//! - exit `2`: usage error. The message names the flag and the offending value.
//!
//! `main` prefixes every message with `mytools <subcommand>:` and prints it on stderr.

use std::fmt;
use std::io;

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    /// Exit 1: an input can't be opened or parsed.
    Data,
    /// Exit 2: the command line itself is invalid.
    Usage,
    /// stdout was closed by the reader (`mytools ... | head`). Not an error for the user.
    BrokenPipe,
}

#[derive(Debug)]
pub struct Error {
    kind: Kind,
    message: String,
    usage: Option<&'static str>,
}

impl Error {
    /// A usage error (exit 2).
    pub fn usage(message: impl Into<String>) -> Self {
        Error { kind: Kind::Usage, message: message.into(), usage: None }
    }

    /// A data error without a location (exit 1).
    pub fn data(message: impl Into<String>) -> Self {
        Error { kind: Kind::Data, message: message.into(), usage: None }
    }

    /// A data error at a 1-based line of a file (exit 1): `file:line: message`.
    pub fn at(file: &str, line: u64, message: impl fmt::Display) -> Self {
        Error::data(format!("{file}:{line}: {message}"))
    }

    /// An input that can't be opened (exit 1): `cannot open file: reason`.
    pub fn open(file: &str, err: &io::Error) -> Self {
        Error::data(format!("cannot open {file}: {}", describe_io(err)))
    }

    /// A read error part-way through an input (exit 1).
    pub fn read(file: &str, err: &io::Error) -> Self {
        Error::data(format!("{file}: {}", describe_io(err)))
    }

    /// A failure writing to stdout. A closed pipe is reported as [`Kind::BrokenPipe`].
    pub fn write(err: io::Error) -> Self {
        if err.kind() == io::ErrorKind::BrokenPipe {
            Error { kind: Kind::BrokenPipe, message: String::new(), usage: None }
        } else {
            Error::data(format!("cannot write output: {}", describe_io(&err)))
        }
    }

    /// Attach the one-line usage string printed after a usage error.
    pub fn with_usage(mut self, usage: &'static str) -> Self {
        self.usage = Some(usage);
        self
    }

    pub fn kind(&self) -> Kind {
        self.kind
    }

    pub fn message(&self) -> &str {
        &self.message
    }

    pub fn usage_line(&self) -> Option<&'static str> {
        self.usage
    }

    /// The process exit code for this error.
    pub fn exit_code(&self) -> u8 {
        match self.kind {
            Kind::Data => 1,
            Kind::Usage => 2,
            Kind::BrokenPipe => 0,
        }
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for Error {}

/// `No such file or directory`, without Rust's ` (os error 2)` suffix.
fn describe_io(err: &io::Error) -> String {
    let text = err.to_string();
    match text.find(" (os error") {
        Some(i) => text[..i].to_string(),
        None => text,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn data_errors_exit_1_and_usage_errors_exit_2() {
        assert_eq!(Error::data("x").exit_code(), 1);
        assert_eq!(Error::usage("x").exit_code(), 2);
    }

    #[test]
    fn located_error_names_file_and_line() {
        let e = Error::at("a.bed", 3, "start 300 is greater than end 200");
        assert_eq!(e.message(), "a.bed:3: start 300 is greater than end 200");
    }

    #[test]
    fn open_error_drops_os_error_suffix() {
        let err = io::Error::from_raw_os_error(2);
        let e = Error::open("nope.bed", &err);
        assert_eq!(e.message(), "cannot open nope.bed: No such file or directory");
        assert_eq!(e.exit_code(), 1);
    }

    #[test]
    fn broken_pipe_is_silent_success() {
        let e = Error::write(io::Error::from(io::ErrorKind::BrokenPipe));
        assert_eq!(e.kind(), Kind::BrokenPipe);
        assert_eq!(e.exit_code(), 0);
    }
}
