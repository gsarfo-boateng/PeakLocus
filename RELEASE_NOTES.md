# PeakLocus 1.0.0 release notes

PeakLocus 1.0.0 is the first packaged release of the RNA-focused peak/interval annotator derived from the original `peak_annotator_modkit_gencode.java` implementation.

The release keeps the original biological goal—mapping RNA-associated genomic coordinates to transcript regions—while generalizing input handling. PeakLocus automatically detects strand columns when possible and safely falls back to unstranded annotation when no reliable strand information exists.

The annotation vocabulary is deliberately limited to RNA-relevant regions:

```text
3'UTR > 5'UTR > CDS > exon > intron > N/A
```

Promoters and other DNA-regulatory annotations are not part of PeakLocus.

This release also derives introns from transcript exon gaps, retains noncoding exons, preserves GENCODE-aware UTR inference, and ships as a zero-dependency Java 8-compatible executable JAR with tests and examples.
