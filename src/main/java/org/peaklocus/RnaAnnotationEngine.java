package org.peaklocus;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;
import java.util.zip.GZIPInputStream;

/**
 * RNA-focused GTF annotation engine used by PeakLocus.
 *
 * Regions are prioritized as:
 *   3'UTR > 5'UTR > CDS > exon > intron > N/A
 *
 * Generic GENCODE "UTR" records are resolved to 5' or 3' UTR from
 * transcript CDS bounds and transcript strand. Introns are inferred from
 * gaps between consecutive exons of each transcript rather than from a
 * generic gene-body overlap.
 */
final class RnaAnnotationEngine {
    private static final Map<String,Integer> PRIORITY = new HashMap<String,Integer>();
    static {
        PRIORITY.put("3'UTR", 1);
        PRIORITY.put("5'UTR", 2);
        PRIORITY.put("CDS", 3);
        PRIORITY.put("exon", 4);
        PRIORITY.put("intron", 5);
        PRIORITY.put("N/A", 6);
    }

    private static final Set<String> INDEX_FEATURES = new HashSet<String>(Arrays.asList(
            "exon", "CDS", "UTR", "start_codon", "stop_codon",
            "five_prime_utr", "three_prime_utr", "five_prime_UTR", "three_prime_UTR",
            "Selenocysteine"));

    static final class CdsBounds {
        int min = Integer.MAX_VALUE;
        int max = Integer.MIN_VALUE;
        boolean seen = false;
        void add(int start, int end) {
            min = Math.min(min, start);
            max = Math.max(max, end);
            seen = true;
        }
    }

    static final class Rec {
        final String chr, strand, feature;
        final String geneId, geneName, transcriptId, geneType, transcriptType;
        final int start, end;

        Rec(String chr, int start, int end, String strand, String feature,
            String geneId, String geneName, String transcriptId,
            String geneType, String transcriptType) {
            this.chr = chr;
            this.start = start;
            this.end = end;
            this.strand = strand;
            this.feature = feature;
            this.geneId = geneId;
            this.geneName = geneName;
            this.transcriptId = transcriptId;
            this.geneType = geneType;
            this.transcriptType = transcriptType;
        }
    }

    static final class Index {
        final int[] starts;
        final int[] prefixMaxEnd;
        final List<Rec> recs;

        Index(List<Rec> input) {
            Collections.sort(input, new Comparator<Rec>() {
                public int compare(Rec a, Rec b) {
                    int c = Integer.compare(a.start, b.start);
                    return c != 0 ? c : Integer.compare(a.end, b.end);
                }
            });
            recs = Collections.unmodifiableList(input);
            starts = new int[input.size()];
            prefixMaxEnd = new int[input.size()];
            int maxEnd = Integer.MIN_VALUE;
            for (int i = 0; i < input.size(); i++) {
                Rec r = input.get(i);
                starts[i] = r.start;
                maxEnd = Math.max(maxEnd, r.end);
                prefixMaxEnd[i] = maxEnd;
            }
        }

        List<Rec> query(int queryStart, int queryEnd) {
            if (starts.length == 0) return Collections.emptyList();

            int lo = 0, hi = starts.length - 1, pos = starts.length;
            while (lo <= hi) {
                int mid = (lo + hi) >>> 1;
                if (starts[mid] > queryEnd) {
                    pos = mid;
                    hi = mid - 1;
                } else {
                    lo = mid + 1;
                }
            }
            if (pos == 0) return Collections.emptyList();

            List<Rec> out = new ArrayList<Rec>();
            for (int i = pos - 1; i >= 0; i--) {
                if (prefixMaxEnd[i] < queryStart) break;
                Rec r = recs.get(i);
                if (r.start <= queryEnd && r.end >= queryStart) out.add(r);
            }
            return out;
        }
    }

    static final class Hit {
        final Rec rec;
        final String location;
        final int overlap;
        final int priority;

        Hit(Rec rec, String location, int overlap) {
            this.rec = rec;
            this.location = location;
            this.overlap = overlap;
            Integer p = PRIORITY.get(location);
            this.priority = p == null ? PRIORITY.get("N/A") : p;
        }
    }

    static final class Result {
        String geneId = "NA";
        String geneName = "NA";
        String transcriptIds = "";
        String geneType = "NA";
        String transcriptType = "NA";
        String geneStrand = "NA";
        String location = "N/A";
    }

    private final Map<String, Map<String, Index>> index;
    private final Map<String, CdsBounds> cdsByTranscript;

    private RnaAnnotationEngine(Map<String, Map<String, Index>> index,
                                Map<String, CdsBounds> cdsByTranscript) {
        this.index = index;
        this.cdsByTranscript = cdsByTranscript;
    }

    static RnaAnnotationEngine load(Path gtf) throws IOException {
        Map<String, CdsBounds> cds = new HashMap<String, CdsBounds>();
        Map<String, Map<String, List<Rec>>> tmp = new HashMap<String, Map<String, List<Rec>>>();
        Map<String, Rec> transcriptMetadata = new HashMap<String, Rec>();
        Map<String, List<Rec>> exonsByTranscript = new HashMap<String, List<Rec>>();

        try (BufferedReader br = open(gtf)) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty() || line.charAt(0) == '#') continue;
                String[] f = line.split("\\t", 9);
                if (f.length < 9) continue;

                String feature = f[2];
                boolean isTranscript = "transcript".equals(feature);
                if (!isTranscript && !INDEX_FEATURES.contains(feature)) continue;

                int start, end;
                try {
                    start = Integer.parseInt(f[3]);
                    end = Integer.parseInt(f[4]);
                } catch (NumberFormatException e) {
                    continue;
                }

                String strand = f[6];
                if (!"+".equals(strand) && !"-".equals(strand)) continue;

                Map<String, String> attr = parseAttrs(f[8]);
                String geneId = value(attr, "gene_id", "");
                String geneName = value(attr, "gene_name", geneId);
                String transcriptId = value(attr, "transcript_id", "");
                String geneType = first(attr, "gene_type", "gene_biotype", "N/A");
                String transcriptType = first(attr, "transcript_type", "transcript_biotype", geneType);

                Rec rec = new Rec(f[0], start, end, strand, feature, geneId, geneName,
                                  transcriptId, geneType, transcriptType);

                if (isTranscript) {
                    if (!transcriptId.isEmpty()) transcriptMetadata.put(transcriptId, rec);
                    continue;
                }

                add(tmp, rec);

                if ("CDS".equals(feature) && !transcriptId.isEmpty()) {
                    CdsBounds b = cds.get(transcriptId);
                    if (b == null) {
                        b = new CdsBounds();
                        cds.put(transcriptId, b);
                    }
                    b.add(start, end);
                }

                if ("exon".equals(feature) && !transcriptId.isEmpty()) {
                    List<Rec> exons = exonsByTranscript.get(transcriptId);
                    if (exons == null) {
                        exons = new ArrayList<Rec>();
                        exonsByTranscript.put(transcriptId, exons);
                    }
                    exons.add(rec);
                    if (!transcriptMetadata.containsKey(transcriptId)) {
                        transcriptMetadata.put(transcriptId,
                                new Rec(rec.chr, rec.start, rec.end, rec.strand, "transcript",
                                        rec.geneId, rec.geneName, rec.transcriptId,
                                        rec.geneType, rec.transcriptType));
                    }
                }
            }
        }

        // Build explicit transcript-aware introns from exon gaps.
        for (Map.Entry<String, List<Rec>> entry : exonsByTranscript.entrySet()) {
            String transcriptId = entry.getKey();
            List<Rec> exons = entry.getValue();
            if (exons.size() < 2) continue;

            Collections.sort(exons, new Comparator<Rec>() {
                public int compare(Rec a, Rec b) {
                    int c = Integer.compare(a.start, b.start);
                    return c != 0 ? c : Integer.compare(a.end, b.end);
                }
            });

            Rec meta = transcriptMetadata.get(transcriptId);
            if (meta == null) meta = exons.get(0);

            int mergedEnd = exons.get(0).end;
            for (int i = 1; i < exons.size(); i++) {
                Rec next = exons.get(i);
                if (next.start > mergedEnd + 1) {
                    Rec intron = new Rec(meta.chr, mergedEnd + 1, next.start - 1,
                                         meta.strand, "intron",
                                         meta.geneId, meta.geneName, transcriptId,
                                         meta.geneType, meta.transcriptType);
                    add(tmp, intron);
                }
                mergedEnd = Math.max(mergedEnd, next.end);
            }
        }

        return new RnaAnnotationEngine(freeze(tmp), cds);
    }

    Result annotate(String chr, int queryStart, int queryEnd, String strand) {
        List<Hit> hits = new ArrayList<Hit>();
        List<String> strands = (strand == null || InputTable.isUnstrandedToken(strand))
                ? Arrays.asList("+", "-")
                : Collections.singletonList(strand);

        for (String s : strands) {
            Map<String, Index> byStrand = index.get(chr);
            if (byStrand == null) continue;
            Index idx = byStrand.get(s);
            if (idx == null) continue;

            for (Rec r : idx.query(queryStart, queryEnd)) {
                int overlap = Math.min(queryEnd, r.end) - Math.max(queryStart, r.start) + 1;
                if (overlap <= 0) continue;
                hits.add(new Hit(r, determineLocation(r), overlap));
            }
        }
        return summarize(hits);
    }

    private Result summarize(List<Hit> hits) {
        Result out = new Result();
        if (hits.isEmpty()) return out;

        Collections.sort(hits, new Comparator<Hit>() {
            public int compare(Hit a, Hit b) {
                int c = Integer.compare(a.priority, b.priority);
                if (c != 0) return c;
                c = Integer.compare(b.overlap, a.overlap);
                if (c != 0) return c;
                c = safe(a.rec.geneId).compareTo(safe(b.rec.geneId));
                if (c != 0) return c;
                return safe(a.rec.transcriptId).compareTo(safe(b.rec.transcriptId));
            }
        });

        Hit best = hits.get(0);
        out.geneId = na(best.rec.geneId);
        out.geneName = na(best.rec.geneName);
        out.geneType = na(best.rec.geneType);
        out.transcriptType = na(best.rec.transcriptType);
        out.geneStrand = na(best.rec.strand);
        out.location = best.location;

        LinkedHashSet<String> transcriptIds = new LinkedHashSet<String>();
        for (Hit h : hits) {
            if (h.priority == best.priority &&
                same(h.location, best.location) &&
                same(h.rec.geneId, best.rec.geneId) &&
                same(h.rec.strand, best.rec.strand) &&
                !safe(h.rec.transcriptId).isEmpty()) {
                transcriptIds.add(h.rec.transcriptId);
            }
        }
        out.transcriptIds = join(transcriptIds, ";");
        return out;
    }

    private String determineLocation(Rec r) {
        String utr = inferUtrType(r);
        if (utr != null) return utr;

        if ("CDS".equals(r.feature) || "Selenocysteine".equals(r.feature)) return "CDS";

        // Kept consistent with the original PeakLocus precursor: codon boundary
        // records are classified with the adjacent UTR category.
        if ("stop_codon".equals(r.feature)) return "3'UTR";
        if ("start_codon".equals(r.feature)) return "5'UTR";

        if ("exon".equals(r.feature)) return "exon";
        if ("intron".equals(r.feature)) return "intron";
        return "N/A";
    }

    private String inferUtrType(Rec r) {
        if ("three_prime_utr".equals(r.feature) || "three_prime_UTR".equals(r.feature)) return "3'UTR";
        if ("five_prime_utr".equals(r.feature) || "five_prime_UTR".equals(r.feature)) return "5'UTR";
        if (!"UTR".equals(r.feature)) return null;

        CdsBounds cds = cdsByTranscript.get(r.transcriptId);
        if (cds == null || !cds.seen) {
            // A generic UTR without CDS bounds cannot be oriented reliably.
            // Treat it as exonic rather than inventing a 5'/3' assignment.
            return "exon";
        }

        if ("+".equals(r.strand)) {
            if (r.end < cds.min) return "5'UTR";
            if (r.start > cds.max) return "3'UTR";
        } else if ("-".equals(r.strand)) {
            if (r.start > cds.max) return "5'UTR";
            if (r.end < cds.min) return "3'UTR";
        }

        // Rare inconsistent annotations: remain conservative.
        return "exon";
    }

    private static void add(Map<String, Map<String, List<Rec>>> map, Rec rec) {
        Map<String, List<Rec>> byStrand = map.get(rec.chr);
        if (byStrand == null) {
            byStrand = new HashMap<String, List<Rec>>();
            map.put(rec.chr, byStrand);
        }
        List<Rec> list = byStrand.get(rec.strand);
        if (list == null) {
            list = new ArrayList<Rec>();
            byStrand.put(rec.strand, list);
        }
        list.add(rec);
    }

    private static Map<String, Map<String, Index>> freeze(Map<String, Map<String, List<Rec>>> tmp) {
        Map<String, Map<String, Index>> out = new HashMap<String, Map<String, Index>>();
        for (Map.Entry<String, Map<String, List<Rec>>> chr : tmp.entrySet()) {
            Map<String, Index> byStrand = new HashMap<String, Index>();
            for (Map.Entry<String, List<Rec>> s : chr.getValue().entrySet()) {
                byStrand.put(s.getKey(), new Index(s.getValue()));
            }
            out.put(chr.getKey(), byStrand);
        }
        return out;
    }

    private static BufferedReader open(Path path) throws IOException {
        InputStream in = Files.newInputStream(path);
        if (path.getFileName().toString().toLowerCase(Locale.ROOT).endsWith(".gz")) {
            in = new GZIPInputStream(in);
        }
        return new BufferedReader(new InputStreamReader(in, StandardCharsets.UTF_8));
    }

    private static Map<String,String> parseAttrs(String attributeString) {
        Map<String,String> out = new HashMap<String,String>();
        for (String attr : attributeString.split(";")) {
            String t = attr.trim();
            if (t.isEmpty()) continue;
            int sp = t.indexOf(' ');
            if (sp > 0) {
                String key = t.substring(0, sp);
                String value = t.substring(sp + 1).replace("\"", "").trim();
                out.put(key, value);
            }
        }
        return out;
    }

    private static String value(Map<String,String> map, String key, String defaultValue) {
        String v = map.get(key);
        return v == null || v.isEmpty() ? defaultValue : v;
    }

    private static String first(Map<String,String> map, String a, String b, String defaultValue) {
        String v = map.get(a);
        if (v != null && !v.isEmpty()) return v;
        v = map.get(b);
        return v == null || v.isEmpty() ? defaultValue : v;
    }

    private static String safe(String s) { return s == null ? "" : s; }
    private static String na(String s) { return s == null || s.isEmpty() ? "NA" : s; }
    private static boolean same(String a, String b) { return a == null ? b == null : a.equals(b); }

    private static String join(Collection<String> values, String sep) {
        StringBuilder sb = new StringBuilder();
        for (String value : values) {
            if (sb.length() > 0) sb.append(sep);
            sb.append(value);
        }
        return sb.toString();
    }
}
