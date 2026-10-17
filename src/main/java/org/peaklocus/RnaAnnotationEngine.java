package org.peaklocus;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;
import java.util.zip.GZIPInputStream;

/**
 * RNA-focused GTF annotation engine used by PeakLocus.
 *
 * Transcript choice is resolved before RNA-location choice. Representative
 * transcript tags and transcript biotype are used to prevent a minor NMD or
 * retained-intron isoform from overriding a canonical/primary coding isoform.
 * Within the selected transcript, specific RNA regions are resolved by the
 * greatest covered number of query bases; exon/intron are structural fallbacks.
 */
final class RnaAnnotationEngine {
    private static final Set<String> INDEX_FEATURES = new HashSet<String>(Arrays.asList(
            "exon", "CDS", "UTR", "start_codon", "stop_codon",
            "five_prime_utr", "three_prime_utr", "five_prime_UTR", "three_prime_UTR",
            "Selenocysteine"));

    private static final List<String> SPECIFIC_LOCATIONS = Arrays.asList("CDS", "5'UTR", "3'UTR");

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

    static final class AttrBag {
        final Map<String,String> values = new HashMap<String,String>();
        final LinkedHashSet<String> tags = new LinkedHashSet<String>();
    }

    static final class TranscriptMeta {
        final String chr, strand, geneId, geneName, transcriptId, geneType, transcriptType;
        final Set<String> tags;
        final int supportLevel;

        TranscriptMeta(String chr, String strand, String geneId, String geneName,
                       String transcriptId, String geneType, String transcriptType,
                       Collection<String> tags, int supportLevel) {
            this.chr = chr;
            this.strand = strand;
            this.geneId = geneId;
            this.geneName = geneName;
            this.transcriptId = transcriptId;
            this.geneType = geneType;
            this.transcriptType = transcriptType;
            this.tags = Collections.unmodifiableSet(new LinkedHashSet<String>(tags));
            this.supportLevel = supportLevel;
        }
    }

    static final class Rec {
        final String chr, strand, feature;
        final String geneId, geneName, transcriptId, geneType, transcriptType;
        final int start, end;
        final Set<String> tags;
        final int transcriptSupportLevel;

        Rec(String chr, int start, int end, String strand, String feature,
            String geneId, String geneName, String transcriptId,
            String geneType, String transcriptType,
            Collection<String> tags, int transcriptSupportLevel) {
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
            this.tags = Collections.unmodifiableSet(new LinkedHashSet<String>(tags));
            this.transcriptSupportLevel = transcriptSupportLevel;
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
        final int overlapStart;
        final int overlapEnd;

        Hit(Rec rec, String location, int queryStart, int queryEnd) {
            this.rec = rec;
            this.location = location;
            this.overlapStart = Math.max(queryStart, rec.start);
            this.overlapEnd = Math.min(queryEnd, rec.end);
            this.overlap = Math.max(0, overlapEnd - overlapStart + 1);
        }
    }

    static final class TranscriptChoice {
        final String transcriptId;
        final Rec representative;
        final List<Hit> hits;
        final int tier;
        final int supportLevel;
        final int locationCoverage;
        final int maxOverlap;
        final String location;

        TranscriptChoice(String transcriptId, Rec representative, List<Hit> hits,
                         int tier, int supportLevel, int locationCoverage,
                         int maxOverlap, String location) {
            this.transcriptId = transcriptId;
            this.representative = representative;
            this.hits = hits;
            this.tier = tier;
            this.supportLevel = supportLevel;
            this.locationCoverage = locationCoverage;
            this.maxOverlap = maxOverlap;
            this.location = location;
        }
    }

    static final class GeneChoice {
        final String geneKey;
        final List<Hit> hits;
        final int maxOverlap;
        final TranscriptChoice transcript;

        GeneChoice(String geneKey, List<Hit> hits, int maxOverlap, TranscriptChoice transcript) {
            this.geneKey = geneKey;
            this.hits = hits;
            this.maxOverlap = maxOverlap;
            this.transcript = transcript;
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
        Map<String, TranscriptMeta> transcriptMeta = loadTranscriptMetadata(gtf);
        Map<String, CdsBounds> cds = new HashMap<String, CdsBounds>();
        Map<String, Map<String, List<Rec>>> tmp = new HashMap<String, Map<String, List<Rec>>>();
        Map<String, List<Rec>> exonsByTranscript = new HashMap<String, List<Rec>>();

        try (BufferedReader br = open(gtf)) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty() || line.charAt(0) == '#') continue;
                String[] f = line.split("\\t", 9);
                if (f.length < 9) continue;

                String feature = f[2];
                if (!INDEX_FEATURES.contains(feature)) continue;

                int start, end;
                try {
                    start = Integer.parseInt(f[3]);
                    end = Integer.parseInt(f[4]);
                } catch (NumberFormatException e) {
                    continue;
                }

                String strand = f[6];
                if (!"+".equals(strand) && !"-".equals(strand)) continue;

                AttrBag attr = parseAttrs(f[8]);
                String transcriptId = value(attr.values, "transcript_id", "");
                TranscriptMeta meta = transcriptMeta.get(transcriptId);

                String geneId = value(attr.values, "gene_id", meta == null ? "" : meta.geneId);
                String geneName = value(attr.values, "gene_name", meta == null ? geneId : meta.geneName);
                String geneType = first(attr.values, "gene_type", "gene_biotype",
                        meta == null ? "N/A" : meta.geneType);
                String transcriptType = first(attr.values, "transcript_type", "transcript_biotype",
                        meta == null ? geneType : meta.transcriptType);

                LinkedHashSet<String> tags = new LinkedHashSet<String>();
                if (meta != null) tags.addAll(meta.tags);
                tags.addAll(attr.tags);

                int supportLevel = parseSupportLevel(value(attr.values, "transcript_support_level", ""));
                if (supportLevel == Integer.MAX_VALUE && meta != null) supportLevel = meta.supportLevel;

                Rec rec = new Rec(f[0], start, end, strand, feature, geneId, geneName,
                        transcriptId, geneType, transcriptType, tags, supportLevel);
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
                }
            }
        }

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

            TranscriptMeta meta = transcriptMeta.get(transcriptId);
            Rec fallback = exons.get(0);
            int mergedEnd = fallback.end;

            for (int i = 1; i < exons.size(); i++) {
                Rec next = exons.get(i);
                if (next.start > mergedEnd + 1) {
                    Rec intron;
                    if (meta != null) {
                        intron = new Rec(meta.chr, mergedEnd + 1, next.start - 1,
                                meta.strand, "intron", meta.geneId, meta.geneName,
                                meta.transcriptId, meta.geneType, meta.transcriptType,
                                meta.tags, meta.supportLevel);
                    } else {
                        intron = new Rec(fallback.chr, mergedEnd + 1, next.start - 1,
                                fallback.strand, "intron", fallback.geneId, fallback.geneName,
                                transcriptId, fallback.geneType, fallback.transcriptType,
                                fallback.tags, fallback.transcriptSupportLevel);
                    }
                    add(tmp, intron);
                }
                mergedEnd = Math.max(mergedEnd, next.end);
            }
        }

        return new RnaAnnotationEngine(freeze(tmp), cds);
    }

    private static Map<String, TranscriptMeta> loadTranscriptMetadata(Path gtf) throws IOException {
        Map<String, TranscriptMeta> out = new HashMap<String, TranscriptMeta>();
        try (BufferedReader br = open(gtf)) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty() || line.charAt(0) == '#') continue;
                String[] f = line.split("\\t", 9);
                if (f.length < 9 || !"transcript".equals(f[2])) continue;

                String strand = f[6];
                if (!"+".equals(strand) && !"-".equals(strand)) continue;

                AttrBag attr = parseAttrs(f[8]);
                String transcriptId = value(attr.values, "transcript_id", "");
                if (transcriptId.isEmpty()) continue;

                String geneId = value(attr.values, "gene_id", "");
                String geneName = value(attr.values, "gene_name", geneId);
                String geneType = first(attr.values, "gene_type", "gene_biotype", "N/A");
                String transcriptType = first(attr.values, "transcript_type", "transcript_biotype", geneType);
                int support = parseSupportLevel(value(attr.values, "transcript_support_level", ""));

                out.put(transcriptId, new TranscriptMeta(f[0], strand, geneId, geneName,
                        transcriptId, geneType, transcriptType, attr.tags, support));
            }
        }
        return out;
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
                Hit hit = new Hit(r, determineLocation(r), queryStart, queryEnd);
                if (hit.overlap > 0) hits.add(hit);
            }
        }
        return summarize(hits);
    }

    private Result summarize(List<Hit> hits) {
        Result out = new Result();
        if (hits.isEmpty()) return out;

        Map<String, List<Hit>> byGene = new HashMap<String, List<Hit>>();
        for (Hit hit : hits) {
            String key = geneKey(hit.rec);
            List<Hit> list = byGene.get(key);
            if (list == null) {
                list = new ArrayList<Hit>();
                byGene.put(key, list);
            }
            list.add(hit);
        }

        List<GeneChoice> genes = new ArrayList<GeneChoice>();
        for (Map.Entry<String, List<Hit>> entry : byGene.entrySet()) {
            int maxOverlap = 0;
            for (Hit h : entry.getValue()) maxOverlap = Math.max(maxOverlap, h.overlap);
            TranscriptChoice tx = chooseBestTranscript(entry.getValue());
            if (tx != null) genes.add(new GeneChoice(entry.getKey(), entry.getValue(), maxOverlap, tx));
        }
        if (genes.isEmpty()) return out;

        Collections.sort(genes, new Comparator<GeneChoice>() {
            public int compare(GeneChoice a, GeneChoice b) {
                int c = Integer.compare(b.maxOverlap, a.maxOverlap);
                if (c != 0) return c;
                c = compareTranscriptChoices(a.transcript, b.transcript);
                if (c != 0) return c;
                return a.geneKey.compareTo(b.geneKey);
            }
        });

        GeneChoice gene = genes.get(0);
        TranscriptChoice best = gene.transcript;
        Rec rec = best.representative;

        out.geneId = na(rec.geneId);
        out.geneName = na(rec.geneName);
        out.geneType = na(rec.geneType);
        out.transcriptType = na(rec.transcriptType);
        out.geneStrand = na(rec.strand);
        out.location = best.location;

        List<TranscriptChoice> allTx = buildTranscriptChoices(gene.hits);
        Collections.sort(allTx, new Comparator<TranscriptChoice>() {
            public int compare(TranscriptChoice a, TranscriptChoice b) {
                return safe(a.transcriptId).compareTo(safe(b.transcriptId));
            }
        });

        LinkedHashSet<String> transcriptIds = new LinkedHashSet<String>();
        if (!safe(best.transcriptId).isEmpty()) transcriptIds.add(best.transcriptId);
        for (TranscriptChoice tx : allTx) {
            if (equivalentTranscriptChoice(tx, best) &&
                same(tx.location, best.location) &&
                same(tx.representative.transcriptType, best.representative.transcriptType) &&
                !safe(tx.transcriptId).isEmpty()) {
                transcriptIds.add(tx.transcriptId);
            }
        }
        out.transcriptIds = join(transcriptIds, ";");
        return out;
    }

    private TranscriptChoice chooseBestTranscript(List<Hit> geneHits) {
        List<TranscriptChoice> choices = buildTranscriptChoices(geneHits);
        if (choices.isEmpty()) return null;
        Collections.sort(choices, new Comparator<TranscriptChoice>() {
            public int compare(TranscriptChoice a, TranscriptChoice b) {
                return compareTranscriptChoices(a, b);
            }
        });
        return choices.get(0);
    }

    private List<TranscriptChoice> buildTranscriptChoices(List<Hit> geneHits) {
        Map<String, List<Hit>> byTranscript = new HashMap<String, List<Hit>>();
        for (Hit hit : geneHits) {
            String tx = safe(hit.rec.transcriptId);
            if (tx.isEmpty()) continue;
            List<Hit> list = byTranscript.get(tx);
            if (list == null) {
                list = new ArrayList<Hit>();
                byTranscript.put(tx, list);
            }
            list.add(hit);
        }

        List<TranscriptChoice> out = new ArrayList<TranscriptChoice>();
        for (Map.Entry<String, List<Hit>> entry : byTranscript.entrySet()) {
            List<Hit> txHits = entry.getValue();
            Rec representative = representativeRec(txHits);
            String location = chooseLocation(txHits);
            int tier = transcriptTier(representative);
            int support = representative.transcriptSupportLevel;
            int locationCoverage = coveredBases(txHits, location);
            int maxOverlap = 0;
            for (Hit h : txHits) maxOverlap = Math.max(maxOverlap, h.overlap);

            out.add(new TranscriptChoice(entry.getKey(), representative, txHits,
                    tier, support, locationCoverage, maxOverlap, location));
        }
        return out;
    }

    private static int compareTranscriptChoices(TranscriptChoice a, TranscriptChoice b) {
        int c = Integer.compare(a.tier, b.tier);
        if (c != 0) return c;
        c = Integer.compare(b.locationCoverage, a.locationCoverage);
        if (c != 0) return c;
        c = Integer.compare(b.maxOverlap, a.maxOverlap);
        if (c != 0) return c;
        c = Integer.compare(a.supportLevel, b.supportLevel);
        if (c != 0) return c;
        return safe(a.transcriptId).compareTo(safe(b.transcriptId));
    }

    private static boolean equivalentTranscriptChoice(TranscriptChoice a, TranscriptChoice b) {
        return a.tier == b.tier &&
               a.locationCoverage == b.locationCoverage &&
               a.maxOverlap == b.maxOverlap &&
               a.supportLevel == b.supportLevel;
    }

    private static Rec representativeRec(List<Hit> hits) {
        Rec best = hits.get(0).rec;
        for (Hit h : hits) {
            Rec r = h.rec;
            int c = Integer.compare(transcriptTier(r), transcriptTier(best));
            if (c < 0 || (c == 0 && r.transcriptSupportLevel < best.transcriptSupportLevel)) best = r;
        }
        return best;
    }

    private String chooseLocation(List<Hit> txHits) {
        int maxSpecificCoverage = 0;
        List<String> winners = new ArrayList<String>();

        for (String location : SPECIFIC_LOCATIONS) {
            int coverage = coveredBases(txHits, location);
            if (coverage > maxSpecificCoverage) {
                maxSpecificCoverage = coverage;
                winners.clear();
                winners.add(location);
            } else if (coverage > 0 && coverage == maxSpecificCoverage) {
                winners.add(location);
            }
        }

        if (winners.size() == 1) return winners.get(0);
        if (winners.size() > 1) return "exon";
        if (coveredBases(txHits, "exon") > 0) return "exon";
        if (coveredBases(txHits, "intron") > 0) return "intron";
        return "N/A";
    }

    private static int coveredBases(List<Hit> hits, String location) {
        List<int[]> intervals = new ArrayList<int[]>();
        for (Hit h : hits) {
            if (same(location, h.location) && h.overlap > 0) {
                intervals.add(new int[]{h.overlapStart, h.overlapEnd});
            }
        }
        if (intervals.isEmpty()) return 0;

        Collections.sort(intervals, new Comparator<int[]>() {
            public int compare(int[] a, int[] b) {
                int c = Integer.compare(a[0], b[0]);
                return c != 0 ? c : Integer.compare(a[1], b[1]);
            }
        });

        int total = 0;
        int start = intervals.get(0)[0];
        int end = intervals.get(0)[1];
        for (int i = 1; i < intervals.size(); i++) {
            int[] x = intervals.get(i);
            if (x[0] <= end + 1) {
                end = Math.max(end, x[1]);
            } else {
                total += end - start + 1;
                start = x[0];
                end = x[1];
            }
        }
        total += end - start + 1;
        return total;
    }

    private static int transcriptTier(Rec r) {
        if (hasAnyTag(r.tags, "gencode_primary", "ensembl_canonical", "mane_select", "mane_plus_clinical")) {
            return 0;
        }
        if (hasTagPrefix(r.tags, "appris_principal")) return 1;

        String type = safe(r.transcriptType).toLowerCase(Locale.ROOT);
        if ("protein_coding".equals(type)) return 2;
        if ("nonsense_mediated_decay".equals(type) || "non_stop_decay".equals(type)) return 4;
        if ("retained_intron".equals(type) || "processed_transcript".equals(type) ||
            "protein_coding_cds_not_defined".equals(type)) return 5;
        return 3;
    }

    private static boolean hasAnyTag(Collection<String> tags, String... wanted) {
        Set<String> normalized = new HashSet<String>();
        for (String tag : tags) normalized.add(tag.toLowerCase(Locale.ROOT));
        for (String w : wanted) if (normalized.contains(w.toLowerCase(Locale.ROOT))) return true;
        return false;
    }

    private static boolean hasTagPrefix(Collection<String> tags, String prefix) {
        String p = prefix.toLowerCase(Locale.ROOT);
        for (String tag : tags) {
            if (tag.toLowerCase(Locale.ROOT).startsWith(p)) return true;
        }
        return false;
    }

    private String determineLocation(Rec r) {
        String utr = inferUtrType(r);
        if (utr != null) return utr;

        if ("CDS".equals(r.feature) || "Selenocysteine".equals(r.feature)) return "CDS";

        // Preserve PeakLocus's existing codon-boundary convention.
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
        if (cds == null || !cds.seen) return "exon";

        if ("+".equals(r.strand)) {
            if (r.end < cds.min) return "5'UTR";
            if (r.start > cds.max) return "3'UTR";
        } else if ("-".equals(r.strand)) {
            if (r.start > cds.max) return "5'UTR";
            if (r.end < cds.min) return "3'UTR";
        }
        return "exon";
    }

    private static String geneKey(Rec r) {
        if (!safe(r.geneId).isEmpty()) return r.geneId;
        if (!safe(r.geneName).isEmpty()) return "NAME:" + r.geneName;
        return "TX:" + safe(r.transcriptId);
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

    private static AttrBag parseAttrs(String attributeString) {
        AttrBag out = new AttrBag();
        for (String attr : attributeString.split(";")) {
            String t = attr.trim();
            if (t.isEmpty()) continue;
            int sp = t.indexOf(' ');
            if (sp <= 0) continue;

            String key = t.substring(0, sp);
            String val = t.substring(sp + 1).replace("\"", "").trim();
            if ("tag".equals(key)) {
                if (!val.isEmpty()) out.tags.add(val);
            } else {
                out.values.put(key, val);
            }
        }
        return out;
    }

    private static int parseSupportLevel(String value) {
        if (value == null) return Integer.MAX_VALUE;
        String s = value.trim();
        if (s.isEmpty() || "NA".equalsIgnoreCase(s)) return Integer.MAX_VALUE;
        int i = 0;
        while (i < s.length() && Character.isDigit(s.charAt(i))) i++;
        if (i == 0) return Integer.MAX_VALUE;
        try {
            return Integer.parseInt(s.substring(0, i));
        } catch (NumberFormatException e) {
            return Integer.MAX_VALUE;
        }
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
