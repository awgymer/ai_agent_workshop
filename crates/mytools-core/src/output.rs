//! Writers for stdout: buffered text (`-nobuf` flushes every line) and BAM, compressed or
//! uncompressed (`-ubam`). stdout carries data only (SPEC §5).

use std::io::{self, BufWriter, StdoutLock, Write};

use noodles::bam;
use noodles::bgzf::io::writer::{Builder, CompressionLevel};
use noodles::sam;

use crate::error::{Error, Result};

/// Line-oriented text output. Each line gets a trailing `\n`.
pub struct TextWriter<W: Write = StdoutLock<'static>> {
    out: BufWriter<W>,
    nobuf: bool,
}

impl TextWriter {
    pub fn stdout(nobuf: bool) -> TextWriter {
        TextWriter::new(io::stdout().lock(), nobuf)
    }
}

impl<W: Write> TextWriter<W> {
    pub fn new(inner: W, nobuf: bool) -> TextWriter<W> {
        TextWriter { out: BufWriter::with_capacity(64 * 1024, inner), nobuf }
    }

    /// Write `text` and a newline.
    pub fn line(&mut self, text: &str) -> Result<()> {
        self.out.write_all(text.as_bytes()).map_err(Error::write)?;
        self.end_line()
    }

    /// Write a line built from parts, without allocating: `w.fmt_line(format_args!(..))`.
    pub fn fmt_line(&mut self, args: std::fmt::Arguments) -> Result<()> {
        self.out.write_fmt(args).map_err(Error::write)?;
        self.end_line()
    }

    /// Write part of a line; finish it with [`TextWriter::end_line`].
    pub fn write(&mut self, text: &str) -> Result<()> {
        self.out.write_all(text.as_bytes()).map_err(Error::write)
    }

    pub fn end_line(&mut self) -> Result<()> {
        self.out.write_all(b"\n").map_err(Error::write)?;
        if self.nobuf {
            self.out.flush().map_err(Error::write)?;
        }
        Ok(())
    }

    /// Flush and hand back the underlying writer.
    pub fn finish(self) -> Result<W> {
        self.out.into_inner().map_err(|e| Error::write(e.into_error()))
    }
}

/// BAM output. The header is written on creation.
pub struct BamWriter<W: Write = StdoutLock<'static>> {
    inner: bam::io::Writer<noodles::bgzf::io::Writer<W>>,
    header: sam::Header,
}

impl BamWriter {
    pub fn stdout(header: &sam::Header, uncompressed: bool) -> Result<BamWriter> {
        BamWriter::new(io::stdout().lock(), header, uncompressed)
    }
}

impl<W: Write> BamWriter<W> {
    pub fn new(inner: W, header: &sam::Header, uncompressed: bool) -> Result<BamWriter<W>> {
        let mut builder = Builder::default();
        if uncompressed {
            builder = builder.set_compression_level(CompressionLevel::NONE);
        }
        let mut inner = bam::io::Writer::from(builder.build_from_writer(inner));
        inner.write_header(header).map_err(Error::write)?;
        Ok(BamWriter { inner, header: header.clone() })
    }

    pub fn write(&mut self, record: &bam::Record) -> Result<()> {
        self.inner.write_record(&self.header, record).map_err(Error::write)
    }

    /// Write the BGZF end-of-file block, flush, and hand back the underlying writer.
    pub fn finish(self) -> Result<W> {
        let mut inner = self.inner.into_inner().finish().map_err(Error::write)?;
        inner.flush().map_err(Error::write)?;
        Ok(inner)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lines_end_with_newline() {
        let mut w = TextWriter::new(Vec::new(), false);
        w.line("chr1\t0\t100").unwrap();
        w.write("chr1\t").unwrap();
        w.write("5").unwrap();
        w.end_line().unwrap();
        w.fmt_line(format_args!("{}\t{}", "chr2", 7)).unwrap();
        assert_eq!(w.finish().unwrap(), b"chr1\t0\t100\nchr1\t5\nchr2\t7\n");
    }

    const BGZF_EOF: [u8; 28] = [
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
        0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ];

    #[test]
    fn bam_writer_output_reads_back_as_bam() {
        for uncompressed in [false, true] {
            let buf = BamWriter::new(Vec::new(), &sam::Header::default(), uncompressed).unwrap().finish().unwrap();
            assert!(buf.ends_with(&BGZF_EOF), "uncompressed={uncompressed}");
            let mut input = crate::input::Input::from_reader("t.bam", Box::new(io::Cursor::new(buf)), None).unwrap();
            assert_eq!(input.format(), crate::record::Format::Bam);
            assert!(input.next_record().unwrap().is_none());
        }
    }
}
