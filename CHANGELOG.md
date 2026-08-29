# Changelog

## PeakLocus 1.0.0 — 2026-08-29

Initial PeakLocus release.

### RNA-focused annotation

- Renamed the software to **PeakLocus**.
- Restricted annotation to RNA/transcript features: `3'UTR`, `5'UTR`, `CDS`, `exon`, `intron`, and `N/A`.
- Removed promoter and DNA-regulatory annotation logic.
- Retained GENCODE generic-UTR inference using CDS bounds and transcript strand.
- Added transcript-aware intron derivation from exon gaps.
- Retained noncoding exon annotation.

### Input handling

- Automatic strand-column detection from headers, BED6/MACS conventions, or strand-value profiles.
- Safe unstranded fallback that queries both transcript strands.
- BED/TSV/CSV/whitespace input support.
- BED3-BED12, modkit, MACS narrowPeak/broadPeak, and custom interval support.
- Automatic or explicit coordinate-system handling.
- Plain-text and gzip input support.

### Packaging

- Executable Java 8-compatible JAR.
- Zero third-party dependencies.
- Build and run scripts.
- Synthetic examples and GTF reference.
- Unit-style and integration tests.
- Detailed README.
