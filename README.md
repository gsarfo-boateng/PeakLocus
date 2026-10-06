# PeakLocus

[![Version](https://img.shields.io/badge/version-0.1.0.0-blue.svg)](#citation)
[![Java](https://img.shields.io/badge/Java-8%2B-orange.svg)](#installation)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

**PeakLocus** is a zero-dependency Java command-line tool for **RNA-centric genomic peak and interval annotation**. It maps genomic intervals from RNA-associated assays to transcript annotations and reports gene identity, transcript identity, biotype, strand, and RNA-region context.

Version: **1.0.0**

<img width="1122" height="1402" alt="PeakLocus" src="https://github.com/user-attachments/assets/eb552016-b9ef-4368-9add-5ee165172211" />


## Purpose

PeakLocus was developed for annotating RNA-associated genomic coordinates rather than general DNA regulatory elements. Its annotation vocabulary is intentionally restricted to transcript-relevant regions:

```text
3'UTR
5'UTR
CDS
exon
intron
N/A
```

PeakLocus does **not** annotate promoters, enhancers, nearest TSSs, chromatin states, or other DNA-regulatory features.

Typical inputs include RNA modification calls, direct RNA/Nanopore modification sites, MeRIP/m6A-seq peaks, CLIP/RBP peaks, RNA-associated MACS2/MACS3 peak calls, and custom BED/TSV interval files.

## What PeakLocus adds

Every input row is preserved and the following columns are appended:

- `Gene_ID`
- `Gene_Name`
- `Transcript_IDs`
- `Gene_Type`
- `Transcript_Type`
- `Gene_Strand`
- `Location`

`Location` is assigned with this priority:

```text
3'UTR > 5'UTR > CDS > exon > intron > N/A
```

This priority is useful when one interval overlaps multiple transcript features. For example, an exon that also has a specific 3'UTR annotation is reported as `3'UTR` rather than only `exon`.

## GENCODE-aware RNA annotation

GENCODE commonly uses the generic feature label `UTR` rather than always writing `five_prime_utr` and `three_prime_utr`. PeakLocus uses transcript CDS boundaries together with transcript strand to infer whether a generic UTR is 5'UTR or 3'UTR.

For a `+` transcript:

```text
5'UTR ---- CDS ---------------- 3'UTR
```

For a `-` transcript, genomic orientation is reversed and PeakLocus accounts for this automatically.

Noncoding transcript exons are retained as `exon`. If a generic UTR cannot be oriented because a transcript has no CDS bounds, PeakLocus conservatively reports it as `exon` rather than inventing a 5'/3' assignment.

## Transcript-aware introns

PeakLocus derives introns from the gaps between consecutive exons of each transcript. This is more RNA-specific than labeling every gene-body region as intronic and prevents an explicitly exonic site from being called intronic merely because it lies inside a gene boundary.

## Automatic strand-column detection

The input does **not** need to place strand in a fixed column.

With the default:

```text
--strand-column auto
```

PeakLocus uses the following logic:

1. **Header detection.** Headers containing `strand` are preferred, including names such as `strand`, `genomic_strand`, `peak_strand`, `read_strand`, or `transcript_strand`.
2. **BED6/MACS convention.** If column 6 contains strand-like values, it is recognized as the standard BED/MACS strand column.
3. **Value-profile detection.** For headerless/custom tables, PeakLocus profiles columns and can identify an arbitrary column dominated by `+` and `-` values.
4. **Safe unstranded fallback.** If no reliable strand column exists, PeakLocus does not invent one. It treats the interval as unstranded and queries both `+` and `-` transcript annotations.

Accepted strand values include:

```text
+  -
plus  minus
forward  reverse
fwd  rev
sense  antisense
```

Values such as `.`, `*`, `NA`, `N/A`, `unknown`, and `unstranded` are treated as unstranded for that row.

At startup PeakLocus reports what it detected, for example:

```text
Strand column: 5
Strand detection: header-based detection: 'strand'
```

You can override detection if necessary:

```bash
--strand-column 4
--strand-column genomic_strand
--strand-column none
```

Column numbers supplied on the command line are 1-based.

## Supported input styles

PeakLocus accepts tab-delimited, comma-delimited, or whitespace-delimited genomic interval tables. Plain text and `.gz` inputs are supported.

Examples include:

- BED3/BED4/BED5/BED6/BED12
- modkit BED/TSV output
- MACS2/MACS3 `narrowPeak`
- MACS2/MACS3 `broadPeak`
- CLIP/RBP peak files
- m6A/MeRIP peak files
- Nanopore/direct RNA modification intervals
- custom RNA peak/interval tables

The input must contain chromosome, start, and end coordinates. With headers, PeakLocus detects common coordinate names automatically. Headerless BED-like files default to columns 1, 2, and 3.

Override unusual layouts with:

```bash
--chrom-column NAME_OR_NUMBER
--start-column NAME_OR_NUMBER
--end-column NAME_OR_NUMBER
```

## Coordinate systems

GTF coordinates are 1-based inclusive. BED and MACS coordinates are normally 0-based half-open.

PeakLocus defaults to:

```text
--coordinates auto
```

It recognizes BED-style names such as `chromStart/chromEnd`, BED file extensions, and `narrowPeak`/`broadPeak`. BED intervals are converted internally from `[start,end)` to GTF-compatible `[start+1,end]`; the original input columns are never modified in the output.

For a custom 1-based inclusive input table:

```bash
--coordinates one-based
```

## Requirements

- Java 8 or newer
- A GTF annotation file such as GENCODE
- No external Java dependencies

## Quick start

```bash
java -jar peaklocus-1.0.0.jar -i modification_sites.bed -g gencode.vM39.primary_assembly.annotation.gtf -o modification_sites.peaklocus.tsv
```

Or use the included wrapper:

```bash
./run.sh -i modification_sites.bed -g gencode.vM39.primary_assembly.annotation.gtf
```

## Examples

### modkit-style table with strand in column 5

```bash
java -jar peaklocus-1.0.0.jar -i examples/input/modkit_like.tsv -g examples/reference/mini.gtf -o modkit.peaklocus.tsv
```

The header contains `strand`, so PeakLocus detects column 5 automatically.

### Standard BED6

```bash
java -jar peaklocus-1.0.0.jar -i examples/input/bed6.bed -g examples/reference/mini.gtf
```

The strand column is detected as BED column 6.

### Headerless custom file with strand in column 4

```text
chr1 104 110 + 9.5 p1
chr1 879 890 - 8.2 p2
chr1 1049 1060 + 7.1 p3
```

Run:

```bash
java -jar peaklocus-1.0.0.jar -i examples/input/headerless_strand_col4.txt -g examples/reference/mini.gtf
```

PeakLocus profiles the columns and identifies column 4 as strand.

### MACS2/MACS3 narrowPeak

```bash
java -jar peaklocus-1.0.0.jar -i examples/input/macs3.narrowPeak -g examples/reference/mini.gtf
```

MACS peaks commonly contain `.` in BED column 6. PeakLocus recognizes the conventional strand field but treats `.` as unstranded, so both transcript strands are queried.

This input support is intended for **RNA-associated peaks** called with MACS; PeakLocus remains an RNA annotation tool and does not add promoter or chromatin annotations.

## Command-line options

```text
-i, --input FILE
-g, --gtf FILE
-o, --output FILE
--strand-column auto|none|N|NAME
--chrom-column auto|N|NAME
--start-column auto|N|NAME
--end-column auto|N|NAME
--coordinates auto|bed|one-based
--quiet
--version
-h, --help
```

## Building from source

PeakLocus deliberately has no external dependencies.

```bash
./build.sh
```

The build script:

1. compiles the source with Java 8 bytecode compatibility;
2. creates `peaklocus-1.0.0.jar`;
3. compiles the self-tests;
4. runs the self-tests automatically.

### Optional Maven build

A standard `pom.xml` is included for users who prefer Maven:

```bash
mvn package
```

The dependency-free `build.sh` remains the simplest reproducible build and does not require Maven.

For a more detailed description of the annotation logic, see `docs/ALGORITHM.md`.

## Testing

Run the full test suite:

```bash
./tests/run_tests.sh
```

The tests cover:

- modkit-style strand in column 5
- standard BED6 strand detection
- arbitrary header-based strand columns
- headerless strand detection
- BED3 unstranded fallback
- plus-strand 5'UTR/CDS/3'UTR inference
- minus-strand 5'UTR/3'UTR inference
- transcript-aware intron construction
- noncoding exon annotation
- complete command-line integration output

## Output example

Input:

```text
chrom   chromStart   chromEnd   name   strand
chr1    184          190        a      +
```

Output includes:

```text
Gene_ID  Gene_Name  Transcript_IDs  Gene_Type       Transcript_Type  Gene_Strand  Location
G1       Gene1      T1              protein_coding  protein_coding   +            3'UTR
```

## Design scope

PeakLocus intentionally answers an RNA-centered question:

> **Which annotated RNA/transcript region does this genomic interval overlap?**

It is not intended to replace general-purpose chromatin or regulatory-region annotation software.

## License

MIT License. See `LICENSE`.
