//! Opening inputs (SPEC §2): a file, `-` or `stdin`. Compression (gzip, BGZF) and
//! format (BED, GFF/GTF, VCF, BAM) are detected from content, never from the extension.
//!
//! Text formats are read one line at a time and converted to BED coordinates
//! (**0-based, half-open**); the original columns are kept verbatim in each record:
//!
//! - BED: columns 2 and 3 as-is, strand from column 6.
//! - GFF/GTF: start is column 4 minus 1, end is column 5, strand from column 7.
//! - VCF: start is POS minus 1, end is start plus the length of REF.
//! - BAM: start is the alignment start minus 1, end adds the reference span; the record's
//!   text is its BED6 view (name, MAPQ as score, strand). Unmapped reads are skipped; a
//!   mapped read with no position is an error, as bedtools rejects it.
//!
//! Header lines (`#`, `track`, `browser`) before the first record are kept for `-header`.
//! Later ones, and blank lines, are skipped (*measured* on `sort`). `\r\n` is accepted.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Cursor, Read};

use flate2::read::MultiGzDecoder;
use noodles::sam::alignment::Record as _;
use noodles::{bam, sam};

use crate::error::{Error, Result};
use crate::record::{Format, Record, Strand};

const DEFAULT_BUFFER: usize = 64 * 1024;
/// `-iobuf` is honoured up to this size so a large value can't exhaust memory.
const MAX_BUFFER: usize = 16 * 1024 * 1024;

pub struct Input {
    name: String,
    format: Format,
    compressed: bool,
    headers: Vec<String>,
    source: Source,
}

enum Source {
    Text { reader: Box<dyn BufRead>, line: u64, pending: Option<Record> },
    Bam { reader: bam::io::Reader<Box<dyn BufRead>>, header: sam::Header, count: u64 },
}

impl Input {
    /// Open a path, or stdin for `-` and `stdin`. `iobuf` is the `-iobuf` size.
    pub fn open(path: &str, iobuf: Option<u64>) -> Result<Input> {
        if path == "-" || path == "stdin" {
            return Input::from_reader("stdin", Box::new(io::stdin()), iobuf);
        }
        let file = File::open(path).map_err(|e| Error::open(path, &e))?;
        Input::from_reader(path, Box::new(file), iobuf)
    }

    /// Read from any byte stream; `name` is used in error messages.
    pub fn from_reader(name: &str, raw: Box<dyn Read>, iobuf: Option<u64>) -> Result<Input> {
        let capacity = iobuf.map_or(DEFAULT_BUFFER, |n| (n as usize).clamp(1, MAX_BUFFER));
        let read_err = |e: io::Error| Error::read(name, &e);
        let (magic, raw) = peek(raw, 2).map_err(read_err)?;
        // BGZF is a series of gzip members, so one multi-member decoder handles both.
        let compressed = magic == [0x1f, 0x8b];
        let decoded: Box<dyn Read> = if compressed { Box::new(MultiGzDecoder::new(raw)) } else { raw };
        let (magic, decoded) = peek(decoded, 4).map_err(read_err)?;
        let reader: Box<dyn BufRead> = Box::new(BufReader::with_capacity(capacity, decoded));
        let mut input = Input {
            name: name.to_string(),
            format: Format::Bed,
            compressed,
            headers: Vec::new(),
            source: Source::Text { reader, line: 0, pending: None },
        };
        if magic == b"BAM\x01" {
            let Source::Text { reader, .. } = std::mem::replace(&mut input.source, dummy_source()) else {
                unreachable!()
            };
            let mut reader = bam::io::Reader::from(reader);
            let header = reader.read_header().map_err(read_err)?;
            input.format = Format::Bam;
            input.source = Source::Bam { reader, header, count: 0 };
        } else {
            input.read_headers()?;
        }
        Ok(input)
    }

    /// The name used in messages: the path, or `stdin`.
    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn format(&self) -> Format {
        self.format
    }

    /// Whether the input was gzip or BGZF compressed.
    pub fn is_compressed(&self) -> bool {
        self.compressed
    }

    /// Header lines before the first record, without line endings.
    pub fn headers(&self) -> &[String] {
        &self.headers
    }

    /// The SAM header of a BAM input.
    pub fn sam_header(&self) -> Option<&sam::Header> {
        match &self.source {
            Source::Bam { header, .. } => Some(header),
            Source::Text { .. } => None,
        }
    }

    /// The next record, or `None` at the end of the input.
    pub fn next_record(&mut self) -> Result<Option<Record>> {
        match &mut self.source {
            Source::Text { reader, line, pending } => {
                if let Some(record) = pending.take() {
                    return Ok(Some(record));
                }
                while let Some(text) = next_line(reader.as_mut(), line, &self.name)? {
                    if !is_header(&text) {
                        return parse_text(self.format, &self.name, *line, text).map(Some);
                    }
                }
                Ok(None)
            }
            Source::Bam { reader, header, count } => next_bam(reader, header, count, &self.name),
        }
    }

    /// Collect header lines, detect the text format from the first record, and hold that
    /// record back for [`Input::next_record`].
    fn read_headers(&mut self) -> Result<()> {
        let Source::Text { reader, line, pending } = &mut self.source else { return Ok(()) };
        let mut vcf = false;
        while let Some(text) = next_line(reader.as_mut(), line, &self.name)? {
            if is_header(&text) {
                vcf |= text.starts_with("##fileformat=VCF");
                self.headers.push(text);
                continue;
            }
            self.format = if vcf { Format::Vcf } else { detect(&text) };
            *pending = Some(parse_text(self.format, &self.name, *line, text)?);
            break;
        }
        Ok(())
    }
}

fn dummy_source() -> Source {
    Source::Text { reader: Box::new(io::empty()), line: 0, pending: None }
}

/// Read up to `n` bytes and return them along with a reader that still yields them.
fn peek(mut reader: Box<dyn Read>, n: usize) -> io::Result<(Vec<u8>, Box<dyn Read>)> {
    let mut prefix = Vec::with_capacity(n);
    (&mut reader).take(n as u64).read_to_end(&mut prefix)?;
    Ok((prefix.clone(), Box::new(Cursor::new(prefix).chain(reader))))
}

/// The next non-blank line without its line ending; `line` counts every line read.
fn next_line(reader: &mut dyn BufRead, line: &mut u64, name: &str) -> Result<Option<String>> {
    loop {
        let mut text = String::new();
        let n = reader.read_line(&mut text).map_err(|e| Error::read(name, &e))?;
        if n == 0 {
            return Ok(None);
        }
        *line += 1;
        while text.ends_with('\n') || text.ends_with('\r') {
            text.pop();
        }
        if !text.is_empty() {
            return Ok(Some(text));
        }
    }
}

fn is_header(text: &str) -> bool {
    text.starts_with('#') || text.starts_with("track") || text.starts_with("browser")
}

/// Guess the format of a text record: integer columns 2 and 3 make BED; nine columns with
/// integer columns 4 and 5 make GFF; eight or more with an integer column 2 make VCF.
fn detect(text: &str) -> Format {
    let (cols, n) = columns(text);
    let int = |i: usize| i < n && cols[i].parse::<u64>().is_ok();
    if int(1) && int(2) {
        Format::Bed
    } else if n >= 9 && int(3) && int(4) {
        Format::Gff
    } else if n >= 8 && int(1) {
        Format::Vcf
    } else {
        Format::Bed
    }
}

/// The first nine columns, and how many columns there are in total.
fn columns(text: &str) -> ([&str; 9], usize) {
    let mut cols = [""; 9];
    let mut n = 0;
    for col in text.split('\t') {
        if n < cols.len() {
            cols[n] = col;
        }
        n += 1;
    }
    (cols, n)
}

fn parse_text(format: Format, name: &str, line: u64, text: String) -> Result<Record> {
    let (cols, n) = columns(&text);
    let err = |msg: String| Error::at(name, line, msg);
    let int = |i: usize, what: &str| -> Result<u64> {
        cols[i].parse().map_err(|_| err(format!("{what} \"{}\" is not an integer", cols[i])))
    };
    let (start, end, strand) = match format {
        Format::Bed => {
            if n < 3 {
                return Err(err(format!("expected at least 3 columns, found {n}")));
            }
            let strand = if n >= 6 { Strand::parse(cols[5]) } else { Strand::Unknown };
            (int(1, "start")?, int(2, "end")?, strand)
        }
        Format::Gff => {
            if n < 9 {
                return Err(err(format!("expected 9 columns for GFF, found {n}")));
            }
            let start = int(3, "start")?;
            let start = start.checked_sub(1).ok_or_else(|| err("GFF start 0 is not 1-based".into()))?;
            (start, int(4, "end")?, Strand::parse(cols[6]))
        }
        Format::Vcf => {
            if n < 8 {
                return Err(err(format!("expected at least 8 columns for VCF, found {n}")));
            }
            let pos = int(1, "POS")?;
            let start = pos.checked_sub(1).ok_or_else(|| err("VCF POS 0 is not 1-based".into()))?;
            (start, start + cols[3].len() as u64, Strand::Unknown)
        }
        Format::Bam => unreachable!("BAM is not text"),
    };
    if start > end {
        return Err(err(format!("start {start} is greater than end {end}")));
    }
    Ok(Record::new(text, start, end, strand, line))
}

fn next_bam(
    reader: &mut bam::io::Reader<Box<dyn BufRead>>,
    header: &sam::Header,
    count: &mut u64,
    name: &str,
) -> Result<Option<Record>> {
    let read_err = |e: io::Error| Error::read(name, &e);
    loop {
        let mut bam_record = bam::Record::default();
        if reader.read_record(&mut bam_record).map_err(read_err)? == 0 {
            return Ok(None);
        }
        *count += 1;
        let flags = bam_record.flags();
        if flags.is_unmapped() {
            continue;
        }
        let Some(id) = bam_record.reference_sequence_id() else {
            continue;
        };
        let id = id.map_err(read_err)?;
        let chrom = match header.reference_sequences().get_index(id) {
            Some((chrom, _)) => chrom.to_string(),
            None => return Err(Error::at(name, *count, format!("reference sequence {id} is not in the header"))),
        };
        // bedtools reads this as start -1 and rejects the file ("Invalid record"). bedtobam
        // writes one for a zero-length BED record at 0: record 13 of tests/fixtures/a.bam.
        let Some(pos) = bam_record.alignment_start() else {
            return Err(Error::at(name, *count, format!("mapped alignment on {chrom} has no position")));
        };
        let start = usize::from(pos.map_err(read_err)?) as u64 - 1;
        let span = bam_record.alignment_span().transpose().map_err(read_err)?.unwrap_or(0) as u64;
        let mut read_name = bam_record.name().map_or_else(|| "*".to_string(), |n| n.to_string());
        if flags.is_segmented() {
            read_name.push_str(if flags.is_first_segment() { "/1" } else { "/2" });
        }
        let mapq = bam_record.mapping_quality().map_or(255, u8::from);
        let strand = if flags.is_reverse_complemented() { Strand::Minus } else { Strand::Plus };
        let text = format!("{chrom}\t{start}\t{}\t{read_name}\t{mapq}\t{}", start + span, strand.as_char());
        let mut record = Record::new(text, start, start + span, strand, *count);
        record.bam = Some(Box::new(bam_record));
        return Ok(Some(record));
    }
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use super::*;

    fn open(text: &str) -> Input {
        Input::from_reader("t.bed", Box::new(Cursor::new(text.as_bytes().to_vec())), None).unwrap()
    }

    fn all(input: &mut Input) -> Vec<(String, u64, u64, Strand, u64)> {
        let mut out = Vec::new();
        while let Some(r) = input.next_record().unwrap() {
            out.push((r.chrom().to_string(), r.start, r.end, r.strand, r.line));
        }
        out
    }

    fn first_error(text: &str) -> String {
        let mut input = match Input::from_reader("t.bed", Box::new(Cursor::new(text.as_bytes().to_vec())), None) {
            Ok(input) => input,
            Err(e) => return e.message().to_string(),
        };
        loop {
            match input.next_record() {
                Ok(Some(_)) => continue,
                Ok(None) => panic!("no error in {text:?}"),
                Err(e) => return e.message().to_string(),
            }
        }
    }

    #[test]
    fn bed_records_keep_their_columns() {
        let mut input = open("chr1\t300\t400\ta05\t40\t+\nchr1\t0\t100\ta01\t10\t-\n");
        assert_eq!(input.format(), Format::Bed);
        let r = input.next_record().unwrap().unwrap();
        assert_eq!(r.text(), "chr1\t300\t400\ta05\t40\t+");
        assert_eq!((r.start, r.end, r.strand), (300, 400, Strand::Plus));
        assert_eq!(all(&mut input), [("chr1".into(), 0, 100, Strand::Minus, 2)]);
    }

    #[test]
    fn headers_before_first_record_are_kept_later_ones_skipped() {
        let mut input = open("track name=x\n#comment\nbrowser position chr1\nchr1\t5\t10\n#mid\n\nchr1\t1\t3\r\n");
        assert_eq!(input.headers(), ["track name=x", "#comment", "browser position chr1"]);
        let records = all(&mut input);
        assert_eq!(records, [("chr1".into(), 5, 10, Strand::Unknown, 4), ("chr1".into(), 1, 3, Strand::Unknown, 7)]);
    }

    #[test]
    fn zero_length_and_position_zero_are_legal() {
        let mut input = open("chr2\t0\t0\n");
        assert_eq!(all(&mut input), [("chr2".into(), 0, 0, Strand::Unknown, 1)]);
    }

    #[test]
    fn bed_errors_name_file_and_line() {
        assert_eq!(first_error("chr1\t0\t10\nchr1\tabc\t200\n"), "t.bed:2: start \"abc\" is not an integer");
        assert_eq!(first_error("chr1\t100\n"), "t.bed:1: expected at least 3 columns, found 2");
        assert_eq!(first_error("chr1\t300\t200\n"), "t.bed:1: start 300 is greater than end 200");
        assert_eq!(first_error("chr1\t-5\t200\n"), "t.bed:1: start \"-5\" is not an integer");
        assert_eq!(first_error("chr1 5 10\n"), "t.bed:1: expected at least 3 columns, found 1");
    }

    #[test]
    fn gff_is_converted_to_zero_based() {
        let mut input = open("##gff-version 3\nchr17\tHAVANA\tgene\t7668421\t7687490\t.\t-\t.\tID=TP53\n");
        assert_eq!(input.format(), Format::Gff);
        assert_eq!(all(&mut input), [("chr17".into(), 7668420, 7687490, Strand::Minus, 2)]);
    }

    #[test]
    fn vcf_end_is_start_plus_ref_length() {
        let text = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                    chr7\t54722323\t.\tC\tT\t50\tPASS\t.\n\
                    chr7\t54722400\t.\tTAC\tT\t50\tPASS\t.\n";
        let mut input = open(text);
        assert_eq!(input.format(), Format::Vcf);
        assert_eq!(input.headers().len(), 2);
        let records = all(&mut input);
        assert_eq!((records[0].1, records[0].2), (54722322, 54722323));
        assert_eq!((records[1].1, records[1].2), (54722399, 54722402));
    }

    #[test]
    fn gzip_is_detected_from_content() {
        let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        gz.write_all(b"chr1\t0\t100\n").unwrap();
        let bytes = gz.finish().unwrap();
        let mut input = Input::from_reader("t.bed.gz", Box::new(Cursor::new(bytes)), None).unwrap();
        assert!(input.is_compressed());
        assert_eq!(all(&mut input), [("chr1".into(), 0, 100, Strand::Unknown, 1)]);
    }

    fn repo_file(path: &str) -> String {
        format!("{}/../../{path}", env!("CARGO_MANIFEST_DIR"))
    }

    #[test]
    fn bgzipped_vcf_fixture() {
        let mut input = Input::open(&repo_file("data/hg002.vcf.gz"), None).unwrap();
        assert!(input.is_compressed());
        assert_eq!(input.format(), Format::Vcf);
        assert!(input.headers().iter().any(|h| h.starts_with("#CHROM")));
        assert_eq!(all(&mut input).len(), 2436);
    }

    #[test]
    fn bam_fixture_gives_bed6_view() {
        let mut input = Input::open(&repo_file("tests/fixtures/a.sorted.bam"), Some(1 << 20)).unwrap();
        assert_eq!(input.format(), Format::Bam);
        assert!(input.sam_header().is_some());
        let first = input.next_record().unwrap().unwrap();
        assert_eq!(first.text(), "chr1\t0\t100\ta01\t255\t+");
        assert!(first.bam.is_some());
        assert_eq!(all(&mut input).len(), 18);
    }

    #[test]
    fn bam_mapped_read_without_position_is_rejected() {
        let mut input = Input::open(&repo_file("tests/fixtures/a.bam"), None).unwrap();
        let e = loop {
            match input.next_record() {
                Ok(Some(_)) => continue,
                Ok(None) => panic!("a.bam read without error"),
                Err(e) => break e,
            }
        };
        assert_eq!(e.exit_code(), 1);
        assert!(e.message().ends_with(":13: mapped alignment on chr2 has no position"), "{}", e.message());
    }

    #[test]
    fn missing_file_is_a_data_error() {
        let e = Input::open("nope.bed", None).err().unwrap();
        assert_eq!(e.exit_code(), 1);
        assert_eq!(e.message(), "cannot open nope.bed: No such file or directory");
    }

    #[test]
    fn empty_input_has_no_records() {
        let mut input = open("");
        assert!(input.next_record().unwrap().is_none());
    }
}
