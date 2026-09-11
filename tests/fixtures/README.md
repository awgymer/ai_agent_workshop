# Test fixtures

Small, committed, never generated at test time. BED is **0-based, half-open**.
These complement `data/` with one file per input format (SPEC §8, item 2). The commands
below are how each file was made, run from the repo root with bedtools v2.31.1; re-run
them only if a fixture genuinely has to change.

| File | Format | Made with |
|---|---|---|
| `a.bed3` | BED3 | `cut -f1-3 data/a.bed` |
| `a.bed4` | BED4 | `cut -f1-4 data/a.bed` |
| `a.bed.gz` | gzipped BED6 | `gzip -nc data/a.bed` |
| `b.bed.gz` | gzipped BED6 | `gzip -nc data/b.bed` |
| `a.bam` | BAM, unsorted | `bedtools bedtobam -i data/a.bed -g /data/GRCh38.chrom.sizes` |
| `a.sorted.bam` | BAM, sorted, no zero-length records | `bedtools sort -i data/a.bed \| awk -F'\t' '$2<$3' \| bedtools bedtobam -i - -g /data/GRCh38.chrom.sizes` |
| `a.bed12` | BED12, unsorted | hand-written: blocks, a single-block record, overlapping transcripts |
| `genes.gff3` | GFF3, unsorted | hand-written from TP53, WRAP53 and EGFR coordinates in `data/genes.gtf` |
| `hg002.head.vcf` | plain VCF | `zcat data/hg002.vcf.gz \| awk '/^#/ \|\| n++ < 40'` (full header, first 40 records) |
| `genome.txt` | genome file for `-g` | the `chr1 chr2 chr3 chrX chr7 chr17` rows of `/data/GRCh38.chrom.sizes`, in that order |
| `genome.fai` | FASTA index for `-faidx` | the same rows of `/data/GRCh38.fa.fai`, in that order |

`genome.txt` and `genome.fai` deliberately use a non-lexicographic order (`chrX` before
`chr7`, `chr7` before `chr17`) so that `-g` and `-faidx` visibly change the output. They
exist because CI has no `/data`.

## Why there are two BAM files

`a.bam` is the file SPEC §8 asks for, but `data/a.bed` has zero-length intervals and
bedtools cannot read its own BAM of them back: `bedtools intersect -a a.bam ...` and
`bedtools merge -i a.bam` both exit 1 with `Invalid record ... chr2 -1 1 a12`. That makes
`a.bam` a good error-path case (both tools should exit 1) and a poor happy-path case.
`a.sorted.bam` drops the zero-length records and is sorted, so `merge` and the `-sorted`
modes can use it.

## Checking them

Every text fixture reads cleanly with `bedtools sort -i <fixture>`; the BAM files pass
`samtools quickcheck` and `bedtools bamtobed -i <fixture>`.
