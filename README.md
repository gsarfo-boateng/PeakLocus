# PeakLocus

[![Version](https://img.shields.io/badge/version-1.0.0-blue.svg)](#citation)
[![Java](https://img.shields.io/badge/Java-8%2B-orange.svg)](#java-version)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

```text
                                 ____  _____    _    _  ___     _      ___   ____ _   _ ____
                                |  _ \| ____|  / \  | |/ /    | |    / _ \ / ___| | | / ___|
                                | |_) |  _|   / _ \ | ' /     | |   | | | | |   | | | \___ \
                                |  __/| |___ / ___ \| . \     | |___| |_| | |___| |_| |___) |
                                |_|   |_____/_/   \_\_|\_\    |_____|\___/ \____|\___/|____/

                                                         PeakLocus v1.0.0
                                          Transcript-aware RNA-centric genomic peak
                                                  and interval annotation

                                inputs:
                                BED3 | BED6 | BED12 | narrowPeak | broadPeak
                                modkit | TSV/CSV | custom intervals

                                regions:
                                3'UTR | 5'UTR | CDS | exon | intron | N/A

                                platforms:
                                macOS arm64 native | Linux x86_64 native | Java 8+
```

**PeakLocus** is a zero-dependency command-line tool for **RNA-centric genomic peak and interval annotation**.

It maps genomic intervals from RNA-associated assays to transcript annotations and reports gene identity, transcript identity, gene/transcript biotype, strand, and RNA-region context.

Version: **1.0.0** https://www.biorxiv.org/content/10.64898/2026.10.17.759721v1

---

## Purpose

PeakLocus was developed specifically for annotating **RNA-associated genomic coordinates**, rather than general DNA regulatory elements.

Its annotation vocabulary is intentionally restricted to transcript-relevant regions:

```text
3'UTR
5'UTR
CDS
exon
intron
N/A
```

PeakLocus does **not** annotate promoters, enhancers, nearest TSSs, chromatin states, or other DNA-regulatory features.

Typical inputs include:

- RNA modification calls
- direct RNA/Nanopore modification sites
- MeRIP-seq / m6A-seq peaks
- CLIP/RBP peaks
- RNA-associated MACS2/MACS3 peak calls
- modkit BED/TSV output
- custom BED/TSV/CSV genomic intervals

---

## What PeakLocus adds

Every input row is preserved and PeakLocus appends:

```text
Gene_ID
Gene_Name
Transcript_IDs
Gene_Type
Transcript_Type
Gene_Strand
Location
```

For example:

```text
Gene_ID         ENSMUSG000000...
Gene_Name       Cd72
Transcript_IDs  ENSMUST00000107926.9
Gene_Type       protein_coding
Transcript_Type protein_coding
Gene_Strand     -
Location        CDS
```

---

## Representative transcript selection

Genes frequently contain multiple transcript isoforms, and the same genomic interval can correspond to different RNA regions in different transcripts.

For example:

```text
representative protein-coding transcript -> CDS
minor NMD transcript                     -> 3'UTR
```

PeakLocus therefore does **not** globally prioritize RNA-region classes across all overlapping transcripts.

Instead, PeakLocus first determines the most representative transcript context and then assigns the RNA-region location within that transcript.

When the relevant metadata are available in the supplied GTF, transcript models are prioritized approximately as:

```text
GENCODE Primary / Ensembl canonical / MANE
                    |
                    v
             APPRIS principal
                    |
                    v
       other protein-coding transcripts
                    |
                    v
          other transcript models
                    |
                    v
   nonsense-mediated / non-stop decay
                    |
                    v
 retained-intron / processed transcripts
```

The ranking is applied only when the corresponding metadata exist in the annotation.

This prevents a minor transcript isoform from controlling the final annotation simply because it labels the interval as a different RNA-region class.

---

## RNA-region assignment

After representative transcript context has been selected, PeakLocus determines the RNA region overlapping the query interval.

Possible assignments are:

```text
3'UTR
5'UTR
CDS
exon
intron
N/A
```

Within the selected transcript, PeakLocus evaluates query-base coverage across specific transcript regions.

Specific RNA subregions such as:

```text
3'UTR
5'UTR
CDS
```

are retained rather than being reduced to their broader parent `exon` category.

For example:

```text
CDS overlap  = 100 bp
exon overlap = 100 bp
```

is reported as:

```text
CDS
```

because CDS is the more specific annotation for that transcript.

If the interval cannot be resolved to a specific UTR or CDS subclass but overlaps an annotated exon, PeakLocus reports:

```text
exon
```

If the query falls within a transcript-derived intron and no exonic feature applies:

```text
intron
```

If no compatible transcript-region annotation is found:

```text
N/A
```

---

## GENCODE-aware transcript ranking

GENCODE annotations frequently contain transcript-ranking metadata such as:

```text
tag "GENCODE_Primary";
tag "Ensembl_canonical";
tag "MANE_Select";
tag "appris_principal_1";
```

PeakLocus preserves repeated GTF `tag` attributes so multiple transcript-ranking annotations can be used simultaneously.

This is important because a single transcript record can contain multiple `tag` fields.

---

## GENCODE-aware UTR annotation

GENCODE commonly uses the generic feature label:

```text
UTR
```

rather than always providing separate:

```text
five_prime_utr
three_prime_utr
```

records.

PeakLocus uses transcript-specific CDS boundaries together with transcript strand to determine whether a generic UTR corresponds to the 5' or 3' untranslated region.

For a `+` strand transcript:

```text
5'UTR -------- CDS ---------------- 3'UTR
```

For a `-` strand transcript:

```text
3'UTR -------- CDS ---------------- 5'UTR
```

PeakLocus automatically accounts for this reversed genomic orientation.

If a generic UTR cannot be reliably oriented because the transcript lacks usable CDS boundaries, PeakLocus conservatively reports the region as:

```text
exon
```

rather than inventing a 5'/3' assignment.

---

## Transcript-aware introns

PeakLocus constructs introns from gaps between consecutive exons of each transcript independently.

Conceptually:

```text
Exon 1 -------- Intron -------- Exon 2
```

This differs from simply treating all non-exonic sequence within a gene boundary as intronic.

Transcript-aware intron construction is especially important at alternatively spliced loci, where a genomic region may be intronic in one transcript but exonic in another.

---

## Automatic strand-column detection

The input does **not** need to place strand in a fixed column.

With the default:

```text
--strand-column auto
```

PeakLocus uses the following strategy:

1. **Header detection**  
   Headers containing strand-related terms are preferred, including:

   ```text
   strand
   genomic_strand
   peak_strand
   read_strand
   transcript_strand
   ```

2. **BED6/MACS convention**  
   If column 6 contains strand-like values, it is recognized as the conventional BED/MACS strand field.

3. **Value-profile detection**  
   For headerless or custom tables, PeakLocus profiles columns and can identify an arbitrary column dominated by strand-like values.

4. **Safe unstranded fallback**  
   If no reliable strand column is found, PeakLocus does not invent one. The interval is treated as unstranded and both transcript strands are queried.

Accepted positive-strand values include:

```text
+
plus
forward
fwd
sense
```

Accepted negative-strand values include:

```text
-
minus
reverse
rev
antisense
anti-sense
```

Values such as:

```text
.
*
NA
N/A
none
unknown
unstranded
```

are treated as unstranded for that row.

At startup PeakLocus reports what it detected, for example:

```text
Strand column: 5
Strand detection: header-based detection: 'strand'
```

You can override automatic detection:

```bash
--strand-column 4
--strand-column genomic_strand
--strand-column none
```

Column numbers supplied on the command line are **1-based**.

---

## Supported input styles

PeakLocus accepts:

- tab-delimited files
- comma-delimited files
- whitespace-delimited files
- plain-text files
- gzip-compressed `.gz` files

Supported interval styles include:

```text
BED3
BED4
BED5
BED6
BED12
modkit BED/TSV
MACS2/MACS3 narrowPeak
MACS2/MACS3 broadPeak
CLIP/RBP peak files
m6A/MeRIP peak files
Nanopore/direct RNA modification intervals
custom RNA peak/interval tables
```

The input must contain chromosome, start, and end coordinates.

With headers, PeakLocus detects common coordinate names automatically.

Headerless BED-like files default to:

```text
column 1 = chromosome
column 2 = start
column 3 = end
```

Unusual layouts can be specified explicitly:

```bash
--chrom-column NAME_OR_NUMBER
--start-column NAME_OR_NUMBER
--end-column NAME_OR_NUMBER
```

---

## Coordinate systems

GTF coordinates are:

```text
1-based inclusive
```

BED and MACS coordinates are normally:

```text
0-based half-open
```

PeakLocus defaults to:

```text
--coordinates auto
```

BED-like inputs are recognized from common coordinate names, extensions, and file formats.

For example:

```text
BED: [start, end)
```

is converted internally to:

```text
GTF-compatible: [start + 1, end]
```

The original input coordinates are **not modified in the output**.

For custom 1-based inclusive tables:

```bash
--coordinates one-based
```

For explicitly BED-like inputs:

```bash
--coordinates bed
```

---

## Installation

PeakLocus v1.0.0 is distributed as:

```text
macOS Apple Silicon native executable
Linux x86_64 native executable
Portable Java JAR
```

The native executables do **not** require Java.

### macOS Apple Silicon

Download from the PeakLocus v1.0.0 GitHub Release:

```text
peaklocus-1.0.0-macos-arm64.tar.gz
```

Extract:

```bash
tar -xzf peaklocus-1.0.0-macos-arm64.tar.gz
cd peaklocus-1.0.0-macos-arm64
```

Verify:

```bash
./peaklocus --version
```

Expected:

```text
PeakLocus 1.0.0
```

Run:

```bash
./peaklocus -i peaks.bed -g annotation.gtf -o results.tsv
```

If macOS blocks the downloaded executable because it is not notarized:

```bash
xattr -d com.apple.quarantine peaklocus
```

### Linux x86_64

Download:

```text
peaklocus-1.0.0-linux-x86_64.tar.gz
```

Extract:

```bash
tar -xzf peaklocus-1.0.0-linux-x86_64.tar.gz
cd peaklocus-1.0.0-linux-x86_64
```

Verify:

```bash
./peaklocus --version
```

Run:

```bash
./peaklocus -i peaks.bed -g annotation.gtf -o results.tsv
```

### Install PeakLocus in your PATH

To make the native executable available from any directory:

```bash
mkdir -p "$HOME/bin"
cp peaklocus "$HOME/bin/"
chmod +x "$HOME/bin/peaklocus"
```

On macOS with zsh:

```bash
echo 'export PATH="$HOME/bin:$PATH"' >> ~/.zshrc
source ~/.zshrc
```

On Linux with bash:

```bash
echo 'export PATH="$HOME/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

Then run:

```bash
peaklocus -i peaks.bed -g annotation.gtf -o results.tsv
```

---

## Java version

Users who prefer Java can download:

```text
peaklocus-1.0.0.jar
```

Requirements:

```text
Java 8 or newer
```

Run:

```bash
java -jar peaklocus-1.0.0.jar -i peaks.bed -g annotation.gtf -o results.tsv
```

PeakLocus has no external Java runtime dependencies.

---

## Quick start

With a native executable installed:

```bash
peaklocus -i modification_sites.bed -g gencode.vM39.primary_assembly.annotation.gtf -o modification_sites.peaklocus.tsv
```

Or with the Java JAR:

```bash
java -jar peaklocus-1.0.0.jar -i modification_sites.bed -g gencode.vM39.primary_assembly.annotation.gtf -o modification_sites.peaklocus.tsv
```

---

## Examples

### modkit-style table

Native:

```bash
peaklocus -i examples/input/modkit_like.tsv -g examples/reference/mini.gtf -o modkit.peaklocus.tsv
```

Java:

```bash
java -jar peaklocus-1.0.0.jar -i examples/input/modkit_like.tsv -g examples/reference/mini.gtf -o modkit.peaklocus.tsv
```

If the input header contains `strand`, PeakLocus detects it automatically.

### Standard BED6

```bash
peaklocus -i examples/input/bed6.bed -g examples/reference/mini.gtf
```

BED column 6 is automatically recognized as the strand field.

### Headerless file with strand in column 4

Example input:

```text
chr1 104 110 + 9.5 p1
chr1 879 890 - 8.2 p2
chr1 1049 1060 + 7.1 p3
```

PeakLocus can identify column 4 as the strand column by profiling the input.

### MACS2/MACS3 narrowPeak

```bash
peaklocus -i examples/input/macs3.narrowPeak -g examples/reference/mini.gtf
```

MACS peak files commonly contain:

```text
.
```

in BED column 6.

PeakLocus recognizes the conventional strand field but treats `.` as unstranded and queries transcript annotations on both strands.

This support is intended for **RNA-associated peaks** called with MACS. PeakLocus remains an RNA annotation tool and does not add promoter or chromatin annotations.

---

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

Complete usage:

```text
peaklocus -i INPUT -g ANNOTATION.gtf [-o OUTPUT] [options]
```

---

## Building from source

Clone the repository:

```bash
git clone https://github.com/gsarfo-boateng/PeakLocus.git
cd PeakLocus
```

Run the full build:

```bash
./build.sh
```

The build process:

1. compiles the PeakLocus source;
2. maintains Java 8 bytecode compatibility;
3. creates `peaklocus-1.0.0.jar`;
4. compiles the self-tests;
5. runs the test suite automatically;
6. runs the representative-transcript regression test.

### Maven build

```bash
mvn package
```

### Native executable build

Native executable creation requires GraalVM with `native-image`.

```bash
./build-native.sh
```

The native Maven profile produces:

```text
target/peaklocus
```

Native binaries are platform-specific and should be built on the target operating system and architecture.

---

## Testing

Run the full test suite:

```bash
./tests/run_tests.sh
```

The tests cover:

- modkit-style strand detection
- standard BED6 strand detection
- arbitrary header-based strand columns
- headerless strand detection
- BED3 unstranded fallback
- plus-strand UTR/CDS annotation
- minus-strand UTR/CDS annotation
- generic GENCODE UTR interpretation
- transcript-aware intron construction
- noncoding exon annotation
- representative-transcript selection
- canonical protein-coding versus NMD transcript conflicts
- complete command-line integration output

### Representative-transcript regression test

PeakLocus includes a dedicated regression test reproducing:

```text
representative protein-coding transcript -> CDS
NMD transcript -> 3'UTR
```

The expected result is:

```text
Transcript_Type = protein_coding
Location = CDS
```

This test prevents a minor alternative transcript model from overriding the biologically representative protein-coding annotation.

---

## Native/JAR equivalence validation

The macOS ARM64 and Linux x86_64 native executables are built from the same PeakLocus source used for the Java JAR.

For the bundled BED6 test, native and Java outputs are required to be identical.

The Linux GitHub Actions workflow also performs a direct native-versus-JAR comparison before adding the Linux binary to the release.

Example validation:

```text
PASS: JAR and native outputs are byte-for-byte identical
```

---

## Output example

Input:

```text
chrom   chromStart   chromEnd   name   strand
chr1    184          190        a      +
```

PeakLocus appends:

```text
Gene_ID
Gene_Name
Transcript_IDs
Gene_Type
Transcript_Type
Gene_Strand
Location
```

Example:

```text
Gene_ID  Gene_Name  Transcript_IDs  Gene_Type       Transcript_Type  Gene_Strand  Location
G1       Gene1      T1              protein_coding  protein_coding   +            3'UTR
```

---

## Annotation philosophy

PeakLocus intentionally separates three concepts:

```text
Gene identity
Transcript identity
RNA-region location
```

A transcript biotype such as:

```text
protein_coding
nonsense_mediated_decay
lncRNA
snRNA
```

is not itself an RNA-region location.

For example:

```text
Gene_Type = protein_coding
Transcript_Type = protein_coding
Location = CDS
```

and:

```text
Gene_Type = snRNA
Transcript_Type = snRNA
Location = exon
```

are both valid PeakLocus annotations.

---

## Design scope

PeakLocus is designed to answer an RNA-centered question:

> **Which annotated RNA/transcript region does this genomic interval overlap in the most representative transcript context available from the supplied annotation?**

PeakLocus is **not** intended to replace general-purpose chromatin or DNA regulatory annotation software.

It does not attempt to assign:

```text
promoter
enhancer
nearest TSS
chromatin state
intergenic regulatory element
```

unless those concepts are supplied independently by another downstream tool.

---

## Citation

If you use PeakLocus in your work, please cite:

**Boateng-Sarfo G, Pritam S, Signor S. PeakLocus: A transcript-aware framework for RNA-centric genomic annotation. PeakLocus v1.0.0. https://www.biorxiv.org/content/10.64898/2026.10.17.759721v1**

Repository:

```text
https://github.com/gsarfo-boateng/PeakLocus
```

---

## License

PeakLocus is distributed under the **MIT License**.

See [`LICENSE`](LICENSE).
