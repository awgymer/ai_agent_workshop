# SPEC.md — mytools

`mytools` is a Rust reimplementation of five bedtools subcommands. The oracle is real
**bedtools v2.31.1**: for the same arguments and input, `mytools` must produce the same
stdout, except for the deviations listed in §8.

BED is **0-based, half-open**: `chr1 100 200` covers bases 100..199.

Every behaviour marked *measured* below was observed on bedtools v2.31.1 on 2026-09-11.

---

## 1. Scope

**Subcommands in v1:** `sort`, `merge`, `intersect`, `window`, `subtract`.
They share one engine, which is why they were picked:

- record parsing and chromosome ordering (`sort`)
- an overlap test plus a per-chromosome index of `-b` (`intersect`; `window` is
  `intersect` with each `-a` interval widened)
- a one-pass sweep over sorted input (`merge`; `subtract` builds on merged `-b`)

**Full parity.** Every flag listed by `bedtools <cmd> -h` in v2.31.1 (78 flags, §4) and
every input and output format those subcommands accept (§2, §5).

**Explicitly NOT in v1:**

- Every other bedtools subcommand. That includes `closest`, which needs its own
  nearest-neighbour search, tie rules and distance conventions, so it shares the least
  machinery with the rest.
- Flags that bedtools accepts but does not list in `-h`. The `-h` listing is the flag set.
- Any bedtools version other than v2.31.1.
- Copying bedtools' stderr wording. Clearer messages are allowed (§7).
- Copying bedtools when it accepts invalid arguments (§7).

## 2. Input formats

- **Formats:** whatever bedtools v2.31.1 accepts for that argument: BED (BED3 to BED12,
  extra columns carried through), GFF/GTF, VCF and BAM. The usage line in each `-h` is
  the starting list; golden tests decide. Format is detected from file content, not the
  extension, as bedtools does.
- **Libraries:** `noodles` reads and writes BAM, SAM, VCF, GFF and BGZF. A gzip crate
  (e.g. `flate2`) handles plain gzip. Other crates are allowed when needed (§9).
- **File or stdin:** both. `-` and `stdin` both mean stdin (*measured*: `sort -i -` and
  `sort -i stdin` both work).
- **Compressed input:** gzip and BGZF, detected from content. *Measured:* a gzipped `-b`
  gives the same result as the plain file.
- **Header lines** (`track`, `browser`, `#`): skipped by default; printed before results
  with `-header`. *Measured* on `sort`: without `-header`, `track`, `browser` and `#` lines
  are dropped; with `-header`, `track` and `#` lines are printed first.

## 3. Interval semantics

- Coordinate system: **0-based half-open**. This is not a decision; BED says so.
- **Non-zero-length overlap:** `a.start < b.end AND b.start < a.end`, strict `<`.
- **Bookended intervals do not overlap** in `intersect`, `window` or `subtract`.
  *Measured:* `100-200` vs `200-300` gives no hit.
- **Bookended intervals do merge** under `merge -d 0`, the default. *Measured:*
  `100-200` + `200-300` gives `100-300`, while `100-200` + `201-300` stays two records.
  A negative `-d` requires that many bases of overlap.
- **Zero-length intervals** (`start == end`) are legal (*measured*: `sort` accepts
  `chr1 0 0`), and bedtools does **not** apply the predicate above to them. *Measured*
  against A = `100-200` with `intersect -u`: B = `150-150`, `100-100` and `200-200` all
  count as overlapping, including at A's end. Reproduce this; the golden cases define
  the rule, not the predicate.
- **Minimum overlap:** 1 bp by default (`-f`/`-F` default `1E-9`), adjusted by
  `-f`, `-F`, `-r` and `-e` exactly as bedtools does.

## 4. Flags per subcommand

Flag names and meanings match bedtools exactly: single dash, multi-letter (`-wa`,
`-wao`), with the value as a separate argument (`-w 500`). `clap` expects `--long` flags,
so the parser is hand-written (§9).

| Subcommand  | Inputs | Flags in v1 (all of `-h`) | Notes |
|-------------|--------|---------------------------|-------|
| `sort`      | `-i` | `-sizeA` `-sizeD` `-chrThenSizeA` `-chrThenSizeD` `-chrThenScoreA` `-chrThenScoreD` `-g` `-faidx` `-header` (9) | Score sorts need BED5 or wider (bedtools exits 1 on BED3). |
| `merge`     | `-i` | `-s` `-S` `-d` `-c` `-o` `-delim` `-prec` `-bed` `-header` `-nobuf` `-iobuf` (11) | Defaults: `-d 0`, `-c 5`, `-o sum`, `-delim ","`, `-prec 5`. `-S` takes `+` or `-`. |
| `intersect` | `-a`, `-b` | `-wa` `-wb` `-loj` `-wo` `-wao` `-u` `-c` `-C` `-v` `-ubam` `-s` `-S` `-f` `-F` `-r` `-e` `-split` `-g` `-nonamecheck` `-sorted` `-names` `-filenames` `-sortout` `-bed` `-header` `-nobuf` `-iobuf` (27) | `-b` takes several files and wildcards. |
| `window`    | `-a` or `-abam`, `-b` | `-abam` `-ubam` `-bed` `-w` `-l` `-r` `-sw` `-sm` `-Sm` `-u` `-c` `-v` `-header` (13) | `-w`, `-l` and `-r` default to 1000. BAM in gives BAM out unless `-bed`. |
| `subtract`  | `-a`, `-b` | `-A` `-N` `-wb` `-wo` `-s` `-S` `-f` `-F` `-r` `-e` `-split` `-g` `-nonamecheck` `-sorted` `-bed` `-header` `-nobuf` `-iobuf` (18) | |

`merge -o` operations (all 20): `sum` `min` `max` `absmin` `absmax` `mean` `median`
`mode` `antimode` `stdev` `sstdev` `collapse` `distinct` `distinct_sort_num`
`distinct_sort_num_desc` `distinct_only` `count` `count_distinct` `first` `last`.
The rules for pairing several `-c` columns with several `-o` operations follow `merge -h`.

- **Strand-aware flags:** all of them (`-s`, `-S`, `-sm`, `-Sm`, `-sw`).
- **Does `merge` need sorted input?** Yes. It does not sort for you; unsorted input is
  a data error (exit 1, §7). *Measured:* bedtools also rejects it.
- **`-sorted` on `intersect` and `subtract`:** requires input sorted by chrom then start
  and uses the one-pass sweep. `-g` fixes chromosome order across files;
  `-nonamecheck` relaxes chromosome name checks, as in bedtools.
- **`-nobuf`:** flush after every output line. **`-iobuf`:** input buffer size, with
  optional `K`/`M`/`G` suffix. Both are validated, and neither changes stdout content.
- **Does `closest` require sorted input?** Not applicable: `closest` is not in v1.

## 5. Output

- **Text output** (BED, GFF, VCF records, counts): byte-identical to bedtools v2.31.1.
  Tab-separated, `\n` line endings, trailing newline after the last line.
- **Empty result:** nothing on stdout, exit 0.
- **Record order:** as bedtools. `intersect`, `window` and `subtract` keep `-a` order;
  `sort` defines its own.
- **Numbers:** `merge` numeric operations use `-prec` (default 5) and format exactly as
  bedtools does; golden cases pin the formatting.
- **BAM output** (`-abam`, `-ubam`, BAM input without `-bed`): must be the same after
  decoding (same header, same records, same order). Compressed bytes may differ, because
  a different compressor legitimately produces different bytes.
- **`-header`:** print A's header lines before the results.
- stdout carries data only. Every error and warning goes to stderr.

## 6. Memory model and performance

**Targets**, per subcommand, on the same input as bedtools:

| Measure | Target | Kind |
|---|---|---|
| Peak memory (`/usr/bin/time -v`, maximum resident set size) | ≤ bedtools | Hard: a miss is a bug |
| Wall-clock time, median of 5 runs | ≤ bedtools + 5% | Soft: a miss is investigated, it doesn't fail the build |

**Complexity:** no worse than bedtools for each subcommand and each mode where bedtools
changes algorithm. Each subcommand's module documents its chosen algorithm and its time
and memory complexity next to the baseline below.

n = records in `-a`/`-i`, m = records in `-b`, k = output records.

| Subcommand / mode | bedtools baseline | Memory bound |
|---|---|---|
| `sort` | Loads the whole input, then sorts. O(n log n) time. | O(n) |
| `merge` | Requires sorted input; one pass. O(n) time. | O(1) intervals, plus the current group for `collapse`, `distinct`, `median`, `mode` and similar |
| `intersect` (default) | Loads `-b` into a binned per-chromosome index and streams `-a`. Each query costs time proportional to the `-b` records in the bins it touches. | O(m) |
| `intersect -sorted` | One-pass sweep over both sorted inputs. O(n + m + k) time. | O(`-b` records overlapping the current position) |
| `window` | As `intersect` (default), with widened `-a`. No `-sorted` mode. | O(m) |
| `subtract` / `subtract -sorted` | As `intersect` / `intersect -sorted`. | As `intersect` |

These baselines come from bedtools' documentation and design. `bench/run.sh` confirms
them by measurement.

- **Streaming vs in-memory:** `merge` and the `-sorted` modes stream. `sort` and the
  default modes of `intersect`, `window` and `subtract` hold one input in memory, as
  bedtools does.
- **Largest input:** `mytools` sets no limit of its own. Practical limits are `noodles`'
  limits and, for the in-memory modes, RAM, the same as for bedtools. Output is buffered
  unless `-nobuf` is given.
- **Benchmark inputs for v1:** the files already in `data/`: `genes.gtf` (9,364 lines),
  `hg002.vcf.gz` (2,436 records), `hg002.highconf.bed` (387), `genes.bed` (25), `a.bed`
  (22) and `b.bed` (19). `broken.vcf` is excluded: it is invalid on purpose.
- **Limitation:** at these sizes timings measure program start-up and parsing, not the
  algorithm. *Measured:* `bedtools intersect -a data/genes.bed -b data/hg002.vcf.gz -c`
  takes 0.00 s and 9.5 MB. They cannot tell O(n log n) from O(n²). Until a larger fixed
  input exists, performance at scale is estimated from the complexity table above.
- **Reference machine:** the workshop VM (2 cores, 3 GB RAM).

## 7. Errors and exit codes

| Exit | Meaning |
|---|---|
| `0` | Success, including an empty result |
| `1` | Bad input data: the command line was valid, but an input can't be opened or parsed |
| `2` | Usage error: the command line itself is invalid |

bedtools exits 1 for every error, including usage errors; exit 2 is a deliberate
deviation (§8).

- **Messages:** one line on stderr, prefixed `mytools <subcommand>:`. Data errors name
  the file and 1-based line number. Usage errors name the flag and the offending value,
  and say what is allowed. No full help text on error.
- **Usage errors are caught before any input is read**, so stdout is empty.
- **Data errors stop at the first bad record.** Output already written for earlier
  records stays written, because output is streamed.
- **Invalid arguments are rejected even when bedtools accepts them.**

| Situation | bedtools v2.31.1 (*measured*) | mytools stderr (example) | Exit |
|---|---|---|---|
| Success | 0 | — | `0` |
| Malformed BED line (`chr1 abc 200`) | 1; "unable to open file or unable to determine types…" | `a.bed:1: start "abc" is not an integer` | `1` |
| Too few columns | 1; same message | `a.bed:1: expected at least 3 columns, found 2` | `1` |
| `start > end` | 1; `sort`: "Start was greater than end"; `intersect`/`subtract`: type-detection message | `a.bed:1: start 300 is greater than end 200` | `1` |
| Score sort on BED3 | 1 | `a.bed:1: -chrThenScoreA needs BED5 or wider, found 3 columns` | `1` |
| Missing input file | 1; "Unable to open file nope.bed." | `cannot open nope.bed: No such file or directory` | `1` |
| Unsorted input to `merge` or a `-sorted` mode | 1; names the out-of-order record | `a.bed:2: not sorted: chr1:100 comes after chr1:500 (run mytools sort first)` | `1` |
| Unknown flag | 1; full help | `unknown option -zz` + one usage line | `2` |
| No arguments / missing required input | 1; full help | `missing required -a` + one usage line | `2` |
| Unknown subcommand | 1 | `unknown subcommand frobnicate` | `2` |
| Invalid value bedtools rejects: `merge -d abc`, `-S x`, `-prec abc`; `intersect -f abc`, `-f 1.5`, `-f -0.5` | 1; full help | `-f 1.5: must be a number in (0, 1]` | `2` |
| Invalid value bedtools **accepts**: `window -w abc`, `-w -5`, `-w 1.5` | **0**, prints results | `-w abc: must be a non-negative integer` | `2` |
| Unknown `merge -o` operation (`-o bogus`) | **0**, error on stderr, no output | `-o bogus: unknown operation` | `2` |

The last two rows are where "reject invalid input" changes observable behaviour. Add a
row whenever another case turns up.

## 8. Correctness

- **Oracle:** real bedtools **v2.31.1** on the files in `data/`. Non-negotiable. CI
  installs exactly that version and checks `bedtools --version`. If the runner's `apt`
  package isn't 2.31.1, install the pinned release instead.
- **Compared:** stdout and exit code. stderr text is never compared.

**Test plan. v1 is done when items 1–5 pass in CI, the flag-coverage check passes, and
item 6 has been run and its results recorded.**

1. **Golden tests: every flag at least once.** `./tests/run_golden.sh` compares
   `$MYTOOLS <args>` with `bedtools <args>`. A coverage check reads the flag list from
   `bedtools <cmd> -h` and fails if any flag has no case. Add hand-picked combinations
   where the output layout changes: `-wa`/`-wb`, `-wo`/`-wao`/`-loj`, `-s`/`-S` with
   `-u`/`-v`, and `-f`/`-F`/`-r`/`-e`. Commands that need sorted input get it sorted into
   a temporary file first.
2. **Golden tests: every input format at least once.** BED3, BED4, BED6, BED12,
   GFF/GTF, VCF, BAM, gzipped input, and stdin (`-`). New fixtures go in `tests/fixtures/`,
   small and committed, never generated at test time:
   - BAM: `bedtools bedtobam -i data/a.bed -g /data/GRCh38.chrom.sizes`
   - BED12: small, hand-written; confirm bedtools accepts it
   - a small GFF3, a plain (uncompressed) VCF cut from `data/hg002.vcf.gz`, and
     gzipped BED
3. **Unit tests, runnable without bedtools** (`cargo test`). One per edge case:
   bookended, zero-length, nested, position 0, and the overlap test itself. Fixing a
   failing golden test means first adding the unit test that would have caught it.
4. **Tests for the deviations.** Every row in §7 marked exit `2`, plus the data-error
   rows, has a test with a fixed expected exit code and stdout (empty for usage errors).
5. **Randomised comparison with bedtools, fixed seed.** Small random BED inputs with
   coordinates packed into 0–50 on `chr1`/`chr2`, both strands, plenty of bookended,
   nested and zero-length intervals. A few hundred cases per subcommand, compared with
   bedtools. These are Rust tests marked `#[ignore]` (they need bedtools), run with
   `cargo test -- --ignored`. The seed is in the source, so failures reproduce.
6. **Benchmark** (`bench/run.sh`), not a CI check. For each subcommand on the §6 inputs:
   median of 5 runs, wall time and maximum resident set size, `mytools` vs bedtools.
   Memory above bedtools is a bug; time over +5% is reported for investigation.

**Accepted deviations from bedtools:**

| Deviation | Why |
|---|---|
| Usage errors exit `2`, not `1` | Callers can tell a bad command line from bad data. |
| stderr wording differs (one line, names file/line or flag/value) | Clearer messages. |
| Invalid values bedtools accepts (`window -w abc`/`-5`/`1.5`, `merge -o bogus`) are rejected with exit `2` and empty stdout | Quietly treating garbage as valid gives wrong results that look right. |
| BAM output is compared after decoding, not byte for byte | A different compressor legitimately produces different bytes. |

## 9. Language and layout

- **Language:** Rust, for every subcommand and every test except the shell scripts
  `tests/run_golden.sh` and `bench/run.sh`.
- **Dependencies:** `noodles` for genomic formats and a gzip crate (e.g. `flate2`).
  Other crates are allowed when needed; declare them in the workspace `Cargo.toml` with
  a one-line reason.
- **Argument parsing:** hand-written, to match bedtools' single-dash multi-letter flags.

```
Cargo.toml                 # workspace only (no [package]); edition + shared deps
crates/
  mytools-core/            # lib: everything the subcommands share
    src/                   #   record model + noodles readers/writers, chrom ordering,
                           #   overlap test, interval index + sorted sweep, merge,
                           #   error type carrying the exit code
  mytools/                 # bin: builds the `mytools` binary
    src/main.rs            #   argument parsing, dispatch
    src/cmd/{sort,merge,intersect,window,subtract}.rs
    tests/                 #   #[ignore] randomised comparison against bedtools
tests/
  run_golden.sh            # golden tests + flag-coverage check
  fixtures/                # committed BAM, BED12, GFF3, VCF, gzipped fixtures
bench/run.sh               # median of 5, time + max RSS vs bedtools
```

- **Unit tests:** `#[cfg(test)]` modules inside `mytools-core`.
- **Build and run:**
  - `cargo build --release`, then `MYTOOLS=target/release/mytools ./tests/run_golden.sh`
  - or `cargo install --path crates/mytools` to put `mytools` on `PATH`
- **CI order:** Rust toolchain, then pinned bedtools v2.31.1, then
  `cargo build --release`, `cargo test`, `cargo test -- --ignored`, and
  `./tests/run_golden.sh`.
- **`.gitignore`:** add `target/`.

**Where this spec overrides `CLAUDE.md`:** `CLAUDE.md` currently says "no third-party
runtime dependencies, standard library only" and "unit tests live in `tests/`". This
spec allows crates (§2, §9) and puts unit tests in `#[cfg(test)]` modules. Update
`CLAUDE.md` to match.
