# CLAUDE.md

Conventions for the `mytools` project. This **overwrites** the `CLAUDE.md` in the repo
root during setup — the one there is notes for the upstream template, not for your
project. It is read automatically at the start of every session, so edit it as your
design firms up, either directly or with the `#` prefix from inside a session.

## My fork — fill this in first

**My fork is `awgymer/ai_agent_workshop`.**

Anything that **creates** — `gh issue create`, `gh pr create` — passes
`--repo awgymer/ai_agent_workshop`. Never write to
`SACGF/ai_agent_workshop`: that's the shared upstream template and thirty other people
are working from it.

Reading and reviewing someone else's PR is fine when I name their fork explicitly
(`gh pr diff`, `gh pr review --repo <partner>/ai_agent_workshop`). The rule is about
where new things land, not about what you're allowed to look at.

Pass the flag every time. `gh` with a missing or empty `--repo` does not fail — it
silently resolves to the git remote and exits 0, so a forgotten flag looks exactly like
a success.

## What this is

`mytools` is a small reimplementation of a subset of bedtools: `sort`, `merge`,
`intersect`, `window`, `subtract` (`closest` is out of scope for v1 — see `SPEC.md` §1).
Real `bedtools` v2.31.1 is installed and is the oracle — if our output differs from it
on the same input, we are wrong. `SPEC.md` is the detailed contract.

## Language — fill this in at 0:30

**`mytools` is written in Rust.** Every subcommand and every test. Don't
introduce a second language without asking me.

One codebase, one language: several agents work on this in parallel and they will each
pick their own otherwise. Reimplementing a single subcommand elsewhere is a deliberate
stretch goal, not a default.

## Interval semantics — read this before touching overlap logic

BED is **0-based, half-open**. `chr1 100 200` covers bases 100..199. Therefore:

- Two intervals overlap iff `a.start < b.end AND b.start < a.end`. Note strict `<`.
- Bookended intervals (`a.end == b.start`) do **not** overlap. They do merge under
  `merge -d 0`.
- Zero-length intervals (`start == end`) are legal in our fixtures and bedtools
  handles them in ways you will not guess. Do not "fix" them — match the oracle.

Every off-by-one bug in this project lives in that comparison. When a golden test
fails, look there first.

## Testing

- `./tests/run_golden.sh [subcommand...]` diffs subcommands against real bedtools.
  Cases live in `tests/golden/<subcommand>.sh`; `tests/flag_coverage.sh <subcommand>`
  fails if a flag from `bedtools <subcommand> -h` has no case.
- **Run it before every commit.** It takes seconds; there is no excuse.
- Fixtures are `data/a.bed`, `data/b.bed` (edge cases), `data/genes.bed`, and the
  small format fixtures in `tests/fixtures/` (commands that made them are in its README).
  Do not regenerate or "tidy" them — the edge cases are deliberate.
- New subcommand or flag? Add its golden case in the same commit.
- Unit tests are `#[cfg(test)]` modules (mostly in `crates/mytools-core`) and must run
  without bedtools (`cargo test`). One per edge case: bookended, zero-length, nested,
  position 0, and the overlap predicate itself.
- Fixed a failing golden test? Add the unit test that would have caught it first.
- If bedtools does something surprising, the test encodes bedtools' behaviour.
  Add a comment saying why; do not encode what you think it should do.

## Code

- Prefer streaming I/O. Read line by line, write as you go. `sort` and the default
  (non-`-sorted`) modes of `intersect`, `window` and `subtract` may hold one input in
  memory, as bedtools does (`SPEC.md` §6).
- Read from a file argument or stdin (`-` means stdin).
- Errors go to stderr, never stdout — stdout is data and gets piped.
- Exit codes: `0` success, `1` bad input data, `2` usage error.
- Dependencies: `noodles` (BAM/SAM/VCF/GFF/BGZF) and `flate2` (gzip), per `SPEC.md` §9.
  Other well-known crates are fine when genuinely needed; declare them in the workspace
  `Cargo.toml` with a one-line reason.

## Commits and PRs

- Small commits, one logical change each, message referencing the issue: `sort: handle
  unsorted chrom order (#3)`.
- Branch per issue: `feat/3-sort`. Exception: trivial one-liners go straight to main —
  ask me which I want rather than defaulting to a PR.
- Closing an issue means checking the code does what the issue asked, not remembering
  that you wrote it.
- PRs go to **your own fork** — see the fork rule at the top of this file.
