# PeakLocus annotation algorithm

PeakLocus is intentionally RNA-centric. It annotates an input genomic interval against transcript features in a GTF and chooses the most specific RNA-region interpretation.

## 1. Input schema detection

PeakLocus detects delimiter, header presence, chromosome/start/end columns, coordinate semantics, and strand.

Strand detection proceeds in this order:

1. header names containing `strand` or closely related orientation labels;
2. BED6/MACS column 6 when values are strand-like;
3. column-wise value profiling for `+/-` or recognized strand synonyms;
4. unstranded fallback if no reliable column exists.

A missing or unstranded value causes PeakLocus to query both GTF strands rather than inventing strand information.

## 2. Coordinate normalization

GTF coordinates are 1-based inclusive. BED/MACS intervals are normally 0-based half-open.

For BED input:

```text
[start, end)  ->  [start + 1, end]
```

Only the internal query coordinates are converted; the input columns written to the output are unchanged.

## 3. GTF indexing

PeakLocus indexes RNA-relevant features by chromosome and strand. Records are sorted by start coordinate and accompanied by a prefix maximum-end array. During overlap queries this allows the backward search to stop only when no earlier interval can overlap, including when genes/transcripts are long.

## 4. UTR inference

If the GTF provides explicit `five_prime_utr` or `three_prime_utr`, PeakLocus uses that label directly.

For generic GENCODE `UTR` records, transcript CDS bounds are calculated first.

On the `+` strand:

```text
UTR left of CDS  ->  5'UTR
UTR right of CDS ->  3'UTR
```

On the `-` strand the interpretation is reversed.

If a generic UTR has no usable CDS bounds, PeakLocus conservatively labels it `exon`.

## 5. Intron construction

Introns are derived per transcript from gaps between consecutive exon intervals. This avoids using the entire gene body as a proxy for intronic sequence and keeps the annotation transcript-aware.

## 6. Region priority

When an interval overlaps more than one RNA feature, PeakLocus chooses the highest-priority category:

```text
3'UTR > 5'UTR > CDS > exon > intron > N/A
```

Within the same category, the hit with the largest base-pair overlap is preferred. Transcript IDs belonging to the selected gene, strand, and location are retained together.

## 7. Output

The original row is preserved and PeakLocus appends:

```text
Gene_ID
Gene_Name
Transcript_IDs
Gene_Type
Transcript_Type
Gene_Strand
Location
```
