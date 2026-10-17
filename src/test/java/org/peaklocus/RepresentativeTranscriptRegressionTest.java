package org.peaklocus;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;

/** Regression test for canonical/primary transcript selection over an NMD UTR isoform. */
public final class RepresentativeTranscriptRegressionTest {
    public static void main(String[] args) throws Exception {
        Path gtf = Files.createTempFile("peaklocus-representative-transcript", ".gtf");
        try {
            Files.write(gtf, Arrays.asList(
                "chr4\ttest\ttranscript\t100\t500\t.\t-\t.\tgene_id \"G1\"; gene_name \"Cd72\"; gene_type \"protein_coding\"; transcript_id \"CANONICAL_TX\"; transcript_type \"protein_coding\"; tag \"basic\"; tag \"Ensembl_canonical\"; tag \"GENCODE_Primary\"; tag \"appris_principal_1\"; transcript_support_level \"1\";",
                "chr4\ttest\texon\t100\t500\t.\t-\t.\tgene_id \"G1\"; gene_name \"Cd72\"; gene_type \"protein_coding\"; transcript_id \"CANONICAL_TX\"; transcript_type \"protein_coding\";",
                "chr4\ttest\tCDS\t200\t450\t.\t-\t0\tgene_id \"G1\"; gene_name \"Cd72\"; gene_type \"protein_coding\"; transcript_id \"CANONICAL_TX\"; transcript_type \"protein_coding\";",
                "chr4\ttest\ttranscript\t100\t500\t.\t-\t.\tgene_id \"G1\"; gene_name \"Cd72\"; gene_type \"protein_coding\"; transcript_id \"NMD_TX\"; transcript_type \"nonsense_mediated_decay\";",
                "chr4\ttest\texon\t100\t500\t.\t-\t.\tgene_id \"G1\"; gene_name \"Cd72\"; gene_type \"protein_coding\"; transcript_id \"NMD_TX\"; transcript_type \"nonsense_mediated_decay\";",
                "chr4\ttest\tCDS\t400\t450\t.\t-\t0\tgene_id \"G1\"; gene_name \"Cd72\"; gene_type \"protein_coding\"; transcript_id \"NMD_TX\"; transcript_type \"nonsense_mediated_decay\";",
                "chr4\ttest\tUTR\t200\t350\t.\t-\t.\tgene_id \"G1\"; gene_name \"Cd72\"; gene_type \"protein_coding\"; transcript_id \"NMD_TX\"; transcript_type \"nonsense_mediated_decay\";"
            ), StandardCharsets.UTF_8);

            RnaAnnotationEngine engine = RnaAnnotationEngine.load(gtf);
            RnaAnnotationEngine.Result result = engine.annotate("chr4", 250, 350, "-");

            assertEquals("Cd72", result.geneName, "geneName");
            assertEquals("CANONICAL_TX", result.transcriptIds, "transcriptIds");
            assertEquals("protein_coding", result.transcriptType, "transcriptType");
            assertEquals("CDS", result.location, "location");

            System.out.println("Representative transcript regression test: PASS");
            System.out.println("Gene_Name=" + result.geneName);
            System.out.println("Transcript_IDs=" + result.transcriptIds);
            System.out.println("Transcript_Type=" + result.transcriptType);
            System.out.println("Location=" + result.location);
        } finally {
            Files.deleteIfExists(gtf);
        }
    }

    private static void assertEquals(String expected, String actual, String field) {
        if (!expected.equals(actual)) {
            throw new AssertionError(field + ": expected=" + expected + " actual=" + actual);
        }
    }
}
