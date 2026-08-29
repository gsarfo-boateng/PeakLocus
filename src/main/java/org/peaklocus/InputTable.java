package org.peaklocus;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;
import java.util.zip.GZIPInputStream;

/**
 * Detects the structure of BED/TSV/CSV-like genomic interval tables.
 * The most important feature is robust strand-column autodetection.
 */
final class InputTable {
    enum Delimiter { TAB, COMMA, WHITESPACE }
    enum CoordinateMode { AUTO, BED, ONE_BASED }

    static final class Config {
        String chromosomeColumn = "auto";
        String startColumn = "auto";
        String endColumn = "auto";
        String strandColumn = "auto";
        CoordinateMode coordinates = CoordinateMode.AUTO;
        int detectionRows = 5000;
    }

    static final class Schema {
        final Delimiter delimiter;
        final boolean hasHeader;
        final String[] header;
        final int chromosomeIndex;
        final int startIndex;
        final int endIndex;
        final int strandIndex; // -1 means no strand column; intervals are unstranded.
        final CoordinateMode coordinates;
        final String strandDetectionReason;

        Schema(Delimiter delimiter, boolean hasHeader, String[] header,
               int chromosomeIndex, int startIndex, int endIndex, int strandIndex,
               CoordinateMode coordinates, String strandDetectionReason) {
            this.delimiter = delimiter;
            this.hasHeader = hasHeader;
            this.header = header;
            this.chromosomeIndex = chromosomeIndex;
            this.startIndex = startIndex;
            this.endIndex = endIndex;
            this.strandIndex = strandIndex;
            this.coordinates = coordinates;
            this.strandDetectionReason = strandDetectionReason;
        }
    }

    static BufferedReader open(Path path) throws IOException {
        InputStream in = Files.newInputStream(path);
        if (path.getFileName().toString().toLowerCase(Locale.ROOT).endsWith(".gz")) {
            in = new GZIPInputStream(in);
        }
        return new BufferedReader(new InputStreamReader(in, StandardCharsets.UTF_8));
    }

    static boolean isSkippableLine(String line) {
        if (line == null) return true;
        String t = line.trim();
        return t.isEmpty() || t.startsWith("#") || t.startsWith("track ") || t.startsWith("browser ");
    }

    static Schema detect(Path path, Config cfg) throws IOException {
        List<String> sampleLines = new ArrayList<String>();
        try (BufferedReader br = open(path)) {
            String line;
            while ((line = br.readLine()) != null && sampleLines.size() < cfg.detectionRows + 1) {
                if (!isSkippableLine(line)) sampleLines.add(line);
            }
        }
        if (sampleLines.isEmpty()) throw new IOException("Input file has no data lines: " + path);

        Delimiter delimiter = detectDelimiter(sampleLines.get(0));
        List<String[]> rows = new ArrayList<String[]>();
        for (String line : sampleLines) rows.add(split(line, delimiter));

        boolean hasHeader = detectHeader(rows.get(0));
        String[] header = hasHeader ? rows.get(0) : null;
        List<String[]> dataRows = hasHeader ? rows.subList(1, rows.size()) : rows;
        if (dataRows.isEmpty()) throw new IOException("Input contains a header but no data rows: " + path);

        int maxCols = 0;
        for (String[] r : dataRows) maxCols = Math.max(maxCols, r.length);

        int chr = resolveGenomicColumn(cfg.chromosomeColumn, header, 0,
                new String[]{"chrom", "chromosome", "chr", "seqname", "seqnames", "contig"}, maxCols, "chromosome");
        int start = resolveGenomicColumn(cfg.startColumn, header, 1,
                new String[]{"chromstart", "start", "startpos", "start_position", "begin", "pos", "position"}, maxCols, "start");
        int end = resolveGenomicColumn(cfg.endColumn, header, 2,
                new String[]{"chromend", "end", "endpos", "end_position", "stop"}, maxCols, "end");

        StrandDetection strand = resolveStrandColumn(cfg.strandColumn, header, dataRows, maxCols);
        CoordinateMode coord = resolveCoordinates(cfg.coordinates, header, start, end, path);

        return new Schema(delimiter, hasHeader, header, chr, start, end, strand.index, coord, strand.reason);
    }

    static String[] split(String line, Delimiter delimiter) {
        if (delimiter == Delimiter.TAB) return line.split("\\t", -1);
        if (delimiter == Delimiter.COMMA) return splitCsv(line);
        return line.trim().split("\\s+");
    }

    private static String[] splitCsv(String line) {
        List<String> out = new ArrayList<String>();
        StringBuilder cur = new StringBuilder();
        boolean inQuotes = false;
        for (int i = 0; i < line.length(); i++) {
            char c = line.charAt(i);
            if (c == '"') {
                if (inQuotes && i + 1 < line.length() && line.charAt(i + 1) == '"') {
                    cur.append('"'); i++;
                } else inQuotes = !inQuotes;
            } else if (c == ',' && !inQuotes) {
                out.add(cur.toString()); cur.setLength(0);
            } else cur.append(c);
        }
        out.add(cur.toString());
        return out.toArray(new String[out.size()]);
    }

    static String join(String[] fields, Delimiter delimiter) {
        String sep = delimiter == Delimiter.COMMA ? "," : "\t";
        // Normalize whitespace-delimited input to TSV output so appended annotations remain safe.
        if (delimiter == Delimiter.WHITESPACE) sep = "\t";
        if (delimiter != Delimiter.COMMA) return joinSimple(fields, sep);
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < fields.length; i++) {
            if (i > 0) sb.append(',');
            String v = fields[i] == null ? "" : fields[i];
            if (v.indexOf(',') >= 0 || v.indexOf('"') >= 0 || v.indexOf('\n') >= 0) {
                sb.append('"').append(v.replace("\"", "\"\"")).append('"');
            } else sb.append(v);
        }
        return sb.toString();
    }

    private static String joinSimple(String[] fields, String sep) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < fields.length; i++) {
            if (i > 0) sb.append(sep);
            sb.append(fields[i] == null ? "" : fields[i]);
        }
        return sb.toString();
    }

    private static Delimiter detectDelimiter(String line) {
        int tabs = count(line, '\t');
        int commas = count(line, ',');
        if (tabs > 0 && tabs >= commas) return Delimiter.TAB;
        if (commas > 0) return Delimiter.COMMA;
        return Delimiter.WHITESPACE;
    }

    private static int count(String s, char c) {
        int n = 0; for (int i = 0; i < s.length(); i++) if (s.charAt(i) == c) n++; return n;
    }

    private static boolean detectHeader(String[] first) {
        if (first.length < 3) return true;
        // Coordinate columns in BED-like data are numeric in columns 2 and 3.
        if (isInteger(first[1]) && isInteger(first[2])) return false;
        // Strong header tokens override unusual numeric-looking layouts.
        for (String f : first) {
            String n = norm(f);
            if (n.equals("chrom") || n.equals("chromosome") || n.equals("chr") ||
                n.equals("chromstart") || n.equals("chromend") || n.equals("strand") ||
                n.equals("start") || n.equals("end")) return true;
        }
        return true;
    }

    private static int resolveGenomicColumn(String spec, String[] header, int fallback,
                                             String[] candidateNames, int maxCols, String label) {
        if (spec != null && !"auto".equalsIgnoreCase(spec)) {
            return resolveColumnSpec(spec, header, maxCols, label);
        }
        if (header != null) {
            int idx = findHeader(header, candidateNames);
            if (idx >= 0) return idx;
        }
        if (fallback >= maxCols) throw new IllegalArgumentException("Cannot infer " + label + " column: input has only " + maxCols + " columns");
        return fallback;
    }

    private static final class StrandDetection {
        final int index; final String reason;
        StrandDetection(int index, String reason) { this.index = index; this.reason = reason; }
    }

    private static StrandDetection resolveStrandColumn(String spec, String[] header,
                                                        List<String[]> dataRows, int maxCols) {
        if (spec != null && !"auto".equalsIgnoreCase(spec)) {
            if ("none".equalsIgnoreCase(spec) || "unstranded".equalsIgnoreCase(spec)) {
                return new StrandDetection(-1, "forced unstranded by --strand-column " + spec);
            }
            int idx = resolveColumnSpec(spec, header, maxCols, "strand");
            return new StrandDetection(idx, "explicit --strand-column " + spec);
        }

        // 1. Strong header names. Header name must either contain 'strand' or have strand-like data.
        if (header != null) {
            int best = -1;
            double bestScore = -1.0;
            for (int c = 0; c < header.length; c++) {
                String h = norm(header[c]);
                boolean strong = h.equals("strand") || h.contains("strand") || h.equals("strnd");
                boolean weak = h.equals("orientation") || h.equals("sense");
                if (!strong && !weak) continue;
                StrandProfile p = profile(dataRows, c);
                double score = p.recognizedFraction + (strong ? 2.0 : 0.0) + (p.informative > 0 ? 0.5 : 0.0);
                if (score > bestScore) { best = c; bestScore = score; }
            }
            if (best >= 0) return new StrandDetection(best, "header-based detection: '" + header[best] + "'");
        }

        // 2. Standard BED6/narrowPeak/broadPeak convention: column 6 is strand when values are strand-like.
        if (maxCols >= 6) {
            StrandProfile p = profile(dataRows, 5);
            if (p.observed > 0 && p.recognizedFraction >= 0.95 && (p.informative > 0 || p.unstranded > 0)) {
                return new StrandDetection(5, "BED6/narrowPeak convention with strand-like values in column 6");
            }
        }

        // 3. Value-profile detection for arbitrary column order. Require +/- evidence to avoid confusing
        // score/flag columns composed solely of '.', NA, or '*'.
        int best = -1;
        double bestScore = 0.0;
        for (int c = 0; c < maxCols; c++) {
            StrandProfile p = profile(dataRows, c);
            if (p.observed < Math.min(3, dataRows.size())) continue;
            if (p.informative == 0) continue;
            if (p.recognizedFraction < 0.95) continue;
            double informativeFraction = p.observed == 0 ? 0.0 : ((double)p.informative / p.observed);
            double score = p.recognizedFraction + informativeFraction;
            if (score > bestScore) { best = c; bestScore = score; }
        }
        if (best >= 0) {
            return new StrandDetection(best, String.format(Locale.ROOT,
                    "value-profile detection in column %d (strand-like score %.3f)", best + 1, bestScore));
        }

        return new StrandDetection(-1, "no reliable strand column detected; treating intervals as unstranded");
    }

    private static final class StrandProfile {
        int observed, recognized, informative, unstranded;
        double recognizedFraction;
    }

    private static StrandProfile profile(List<String[]> rows, int c) {
        StrandProfile p = new StrandProfile();
        for (String[] r : rows) {
            if (c >= r.length) continue;
            String v = stripQuotes(r[c]).trim();
            if (v.isEmpty()) continue;
            p.observed++;
            String normalized = normalizeStrandToken(v);
            if ("+".equals(normalized) || "-".equals(normalized)) { p.recognized++; p.informative++; }
            else if (normalized == null && isUnstrandedToken(v)) { p.recognized++; p.unstranded++; }
        }
        p.recognizedFraction = p.observed == 0 ? 0.0 : ((double)p.recognized / p.observed);
        return p;
    }


    static String normalizeStrandToken(String v) {
        if (v == null) return null;
        String x = v.trim().toLowerCase(Locale.ROOT);
        if (x.equals("+") || x.equals("plus") || x.equals("forward") || x.equals("fwd") || x.equals("sense")) return "+";
        if (x.equals("-") || x.equals("minus") || x.equals("reverse") || x.equals("rev") || x.equals("antisense") || x.equals("anti-sense")) return "-";
        return null;
    }

    static boolean isUnstrandedToken(String v) {
        if (v == null) return true;
        String x = v.trim().toLowerCase(Locale.ROOT);
        return x.isEmpty() || x.equals(".") || x.equals("*") || x.equals("na") || x.equals("n/a") ||
               x.equals("none") || x.equals("unknown") || x.equals("unstranded");
    }

    private static CoordinateMode resolveCoordinates(CoordinateMode requested, String[] header,
                                                     int startIndex, int endIndex, Path path) {
        if (requested != CoordinateMode.AUTO) return requested;
        if (header != null) {
            String s = startIndex < header.length ? norm(header[startIndex]) : "";
            String e = endIndex < header.length ? norm(header[endIndex]) : "";
            if (s.equals("chromstart") || e.equals("chromend") || s.contains("bedstart")) return CoordinateMode.BED;
            if (s.equals("position") || s.equals("pos") || s.contains("1based")) return CoordinateMode.ONE_BASED;
        }
        String name = path.getFileName().toString().toLowerCase(Locale.ROOT);
        if (name.endsWith(".bed") || name.endsWith(".bed.gz") || name.contains("narrowpeak") || name.contains("broadpeak")) {
            return CoordinateMode.BED;
        }
        // BED semantics are the safest default for RNA peak/interval software.
        return CoordinateMode.BED;
    }

    private static int resolveColumnSpec(String spec, String[] header, int maxCols, String label) {
        try {
            int oneBased = Integer.parseInt(spec);
            int idx = oneBased - 1;
            if (idx < 0 || idx >= maxCols) throw new IllegalArgumentException(label + " column " + spec + " is outside input width " + maxCols);
            return idx;
        } catch (NumberFormatException ignored) { }
        if (header == null) throw new IllegalArgumentException("Cannot resolve " + label + " column name '" + spec + "' because input has no header");
        String target = norm(spec);
        for (int i = 0; i < header.length; i++) if (norm(header[i]).equals(target)) return i;
        throw new IllegalArgumentException("Column '" + spec + "' not found for " + label);
    }

    private static int findHeader(String[] header, String[] names) {
        Set<String> n = new HashSet<String>(); for (String s : names) n.add(norm(s));
        for (int i = 0; i < header.length; i++) if (n.contains(norm(header[i]))) return i;
        return -1;
    }

    private static String norm(String s) {
        if (s == null) return "";
        return stripQuotes(s).trim().toLowerCase(Locale.ROOT).replaceAll("[^a-z0-9]+", "");
    }

    static String stripQuotes(String s) {
        if (s == null) return "";
        String x = s.trim();
        if (x.length() >= 2 && ((x.startsWith("\"") && x.endsWith("\"")) || (x.startsWith("'") && x.endsWith("'")))) {
            return x.substring(1, x.length() - 1);
        }
        return x;
    }

    static boolean isInteger(String s) {
        try { Long.parseLong(stripQuotes(s)); return true; } catch (Exception e) { return false; }
    }
}
