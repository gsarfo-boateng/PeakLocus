package org.peaklocus;

import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;

/** Zero-dependency PeakLocus self-tests. Run with assertions enabled (-ea). */
public final class PeakLocusSelfTest {
    public static void main(String[] args) throws Exception {
        Path dir = Files.createTempDirectory("peaklocus-test-");
        try {
            testModkitHeaderStrand(dir);
            testBed6Strand(dir);
            testCustomHeaderStrand(dir);
            testHeaderlessStrand(dir);
            testBed3Unstranded(dir);
            testPlusStrandUtrAndCds(dir);
            testMinusStrandUtr(dir);
            testTranscriptAwareIntron(dir);
            testNoncodingExon(dir);
            System.out.println("ALL TESTS PASSED");
        } finally {
            deleteTree(dir);
        }
    }

    private static void testModkitHeaderStrand(Path d) throws Exception {
        Path p = d.resolve("modkit.tsv");
        Files.write(p, Arrays.asList(
                "chrom\tchromStart\tchromEnd\tname\tstrand\tvalid_coverage\tpercent_modified",
                "chr1\t100\t101\ta\t+\t20\t80.0",
                "chr1\t184\t185\ta\t-\t25\t76.0"), StandardCharsets.UTF_8);
        InputTable.Schema s = InputTable.detect(p, new InputTable.Config());
        check(s.strandIndex == 4, "modkit-style strand should be detected in column 5");
    }

    private static void testBed6Strand(Path d) throws Exception {
        Path p = d.resolve("x.bed");
        Files.write(p, Arrays.asList(
                "chr1\t0\t10\tp1\t100\t+",
                "chr1\t20\t30\tp2\t100\t-"), StandardCharsets.UTF_8);
        InputTable.Schema s = InputTable.detect(p, new InputTable.Config());
        check(s.strandIndex == 5, "BED6 strand should be column 6");
        check(s.coordinates == InputTable.CoordinateMode.BED, "BED coordinates expected");
    }

    private static void testCustomHeaderStrand(Path d) throws Exception {
        Path p = d.resolve("custom.tsv");
        Files.write(p, Arrays.asList(
                "chrom\tstart\tend\tsignal\tgenomic_strand\tname",
                "chr1\t1\t5\t9.2\t-\tx",
                "chr1\t6\t9\t1.2\t+\ty"), StandardCharsets.UTF_8);
        InputTable.Schema s = InputTable.detect(p, new InputTable.Config());
        check(s.strandIndex == 4, "header strand should be detected in arbitrary column");
    }

    private static void testHeaderlessStrand(Path d) throws Exception {
        Path p = d.resolve("weird.txt");
        Files.write(p, Arrays.asList(
                "chr1 0 10 + 7.0 foo",
                "chr1 20 30 - 3.0 bar",
                "chr1 40 50 + 8.0 baz"), StandardCharsets.UTF_8);
        InputTable.Schema s = InputTable.detect(p, new InputTable.Config());
        check(s.strandIndex == 3, "headerless +/- column should be detected");
    }

    private static void testBed3Unstranded(Path d) throws Exception {
        Path p = d.resolve("x3.bed");
        Files.write(p, Arrays.asList("chr1\t0\t10", "chr1\t20\t30"), StandardCharsets.UTF_8);
        InputTable.Schema s = InputTable.detect(p, new InputTable.Config());
        check(s.strandIndex == -1, "BED3 should be treated as unstranded");
    }

    private static void testPlusStrandUtrAndCds(Path d) throws Exception {
        RnaAnnotationEngine e = RnaAnnotationEngine.load(writeMiniGtf(d));
        RnaAnnotationEngine.Result utr5 = e.annotate("chr1", 105, 110, "+");
        RnaAnnotationEngine.Result cds = e.annotate("chr1", 130, 140, "+");
        RnaAnnotationEngine.Result utr3 = e.annotate("chr1", 185, 190, "+");
        check("5'UTR".equals(utr5.location), "expected inferred 5'UTR");
        check("CDS".equals(cds.location), "expected CDS");
        check("3'UTR".equals(utr3.location), "expected inferred 3'UTR");
    }

    private static void testMinusStrandUtr(Path d) throws Exception {
        RnaAnnotationEngine e = RnaAnnotationEngine.load(writeMiniGtf(d));
        RnaAnnotationEngine.Result high = e.annotate("chr1", 880, 890, "-");
        RnaAnnotationEngine.Result low = e.annotate("chr1", 710, 720, "-");
        check("5'UTR".equals(high.location), "minus-strand high-coordinate UTR should be 5'UTR");
        check("3'UTR".equals(low.location), "minus-strand low-coordinate UTR should be 3'UTR");
    }

    private static void testTranscriptAwareIntron(Path d) throws Exception {
        RnaAnnotationEngine e = RnaAnnotationEngine.load(writeMiniGtf(d));
        RnaAnnotationEngine.Result r = e.annotate("chr1", 250, 260, "+");
        check("intron".equals(r.location), "expected transcript-derived intron");
        check("T1".equals(r.transcriptIds), "expected T1 intron transcript");
    }

    private static void testNoncodingExon(Path d) throws Exception {
        RnaAnnotationEngine e = RnaAnnotationEngine.load(writeMiniGtf(d));
        RnaAnnotationEngine.Result r = e.annotate("chr1", 1050, 1060, "+");
        check("exon".equals(r.location), "noncoding exon should remain exon");
        check("lncRNA".equals(r.geneType), "expected lncRNA gene type");
    }

    private static Path writeMiniGtf(Path d) throws Exception {
        Path g = d.resolve("mini.gtf");
        if (Files.exists(g)) return g;
        Files.write(g, Arrays.asList(
                // + strand protein-coding transcript with two exons and an intron.
                "chr1\ttest\ttranscript\t101\t500\t.\t+\t.\tgene_id \"G1\"; gene_name \"Gene1\"; transcript_id \"T1\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\";",
                "chr1\ttest\texon\t101\t200\t.\t+\t.\tgene_id \"G1\"; gene_name \"Gene1\"; transcript_id \"T1\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\";",
                "chr1\ttest\tUTR\t101\t120\t.\t+\t.\tgene_id \"G1\"; gene_name \"Gene1\"; transcript_id \"T1\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\";",
                "chr1\ttest\tCDS\t121\t180\t.\t+\t0\tgene_id \"G1\"; gene_name \"Gene1\"; transcript_id \"T1\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\";",
                "chr1\ttest\tUTR\t181\t200\t.\t+\t.\tgene_id \"G1\"; gene_name \"Gene1\"; transcript_id \"T1\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\";",
                "chr1\ttest\texon\t301\t500\t.\t+\t.\tgene_id \"G1\"; gene_name \"Gene1\"; transcript_id \"T1\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\";",
                // - strand protein-coding transcript.
                "chr1\ttest\ttranscript\t701\t900\t.\t-\t.\tgene_id \"G2\"; gene_name \"Gene2\"; transcript_id \"T2\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\";",
                "chr1\ttest\texon\t701\t900\t.\t-\t.\tgene_id \"G2\"; gene_name \"Gene2\"; transcript_id \"T2\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\";",
                "chr1\ttest\tUTR\t701\t730\t.\t-\t.\tgene_id \"G2\"; gene_name \"Gene2\"; transcript_id \"T2\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\";",
                "chr1\ttest\tCDS\t731\t870\t.\t-\t0\tgene_id \"G2\"; gene_name \"Gene2\"; transcript_id \"T2\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\";",
                "chr1\ttest\tUTR\t871\t900\t.\t-\t.\tgene_id \"G2\"; gene_name \"Gene2\"; transcript_id \"T2\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\";",
                // Noncoding transcript.
                "chr1\ttest\ttranscript\t1001\t1100\t.\t+\t.\tgene_id \"G3\"; gene_name \"Lnc1\"; transcript_id \"T3\"; gene_type \"lncRNA\"; transcript_type \"lncRNA\";",
                "chr1\ttest\texon\t1001\t1100\t.\t+\t.\tgene_id \"G3\"; gene_name \"Lnc1\"; transcript_id \"T3\"; gene_type \"lncRNA\"; transcript_type \"lncRNA\";"
        ), StandardCharsets.UTF_8);
        return g;
    }

    private static void check(boolean ok, String message) {
        if (!ok) throw new AssertionError(message);
    }

    private static void deleteTree(Path p) throws Exception {
        if (!Files.exists(p)) return;
        List<Path> paths = new ArrayList<Path>();
        Files.walk(p).forEach(paths::add);
        Collections.sort(paths, Collections.reverseOrder());
        for (Path q : paths) Files.deleteIfExists(q);
    }
}
