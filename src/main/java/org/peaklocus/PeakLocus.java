package org.peaklocus;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;

/** Command-line entry point for PeakLocus. */
public final class PeakLocus {
    private static final String VERSION = "1.0.0";
    private static final String[] ANNOTATION_HEADER = new String[]{
            "Gene_ID", "Gene_Name", "Transcript_IDs", "Gene_Type",
            "Transcript_Type", "Gene_Strand", "Location"
    };

    private static final class Cli {
        Path input, gtf, output;
        String chromosomeColumn = "auto";
        String startColumn = "auto";
        String endColumn = "auto";
        String strandColumn = "auto";
        InputTable.CoordinateMode coordinates = InputTable.CoordinateMode.AUTO;
        boolean quiet = false;
    }

    public static void main(String[] args) {
        try {
            Cli cli = parse(args);
            if (cli == null) return;
            run(cli);
        } catch (Exception e) {
            System.err.println("ERROR: " + e.getMessage());
            if (Boolean.getBoolean("peaklocus.debug")) e.printStackTrace(System.err);
            System.exit(2);
        }
    }

    private static void run(Cli cli) throws Exception {
        InputTable.Config cfg = new InputTable.Config();
        cfg.chromosomeColumn = cli.chromosomeColumn;
        cfg.startColumn = cli.startColumn;
        cfg.endColumn = cli.endColumn;
        cfg.strandColumn = cli.strandColumn;
        cfg.coordinates = cli.coordinates;

        InputTable.Schema schema = InputTable.detect(cli.input, cfg);
        if (cli.output == null) cli.output = defaultOutput(cli.input, schema.delimiter);

        if (!cli.quiet) {
            System.err.println("PeakLocus v" + VERSION);
            System.err.println("RNA-centric genomic peak and interval annotation");
            System.err.println("Input: " + cli.input);
            System.err.println("GTF: " + cli.gtf);
            System.err.println("Detected delimiter: " + schema.delimiter);
            System.err.println("Header: " + (schema.hasHeader ? "yes" : "no"));
            System.err.println("Chrom/start/end columns: " + (schema.chromosomeIndex + 1) + "/" +
                               (schema.startIndex + 1) + "/" + (schema.endIndex + 1));
            System.err.println("Strand column: " + (schema.strandIndex >= 0 ?
                               Integer.toString(schema.strandIndex + 1) : "none"));
            System.err.println("Strand detection: " + schema.strandDetectionReason);
            System.err.println("Coordinate mode: " + schema.coordinates);
            System.err.println("Building RNA annotation index...");
        }

        RnaAnnotationEngine engine = RnaAnnotationEngine.load(cli.gtf);

        long rows = 0, written = 0, invalid = 0;
        try (BufferedReader br = InputTable.open(cli.input);
             BufferedWriter bw = Files.newBufferedWriter(cli.output, StandardCharsets.UTF_8)) {

            boolean headerWritten = false;
            String line;
            while ((line = br.readLine()) != null) {
                if (InputTable.isSkippableLine(line)) continue;
                String[] fields = InputTable.split(line, schema.delimiter);

                if (schema.hasHeader && !headerWritten) {
                    bw.write(InputTable.join(append(fields, ANNOTATION_HEADER), schema.delimiter));
                    bw.newLine();
                    headerWritten = true;
                    continue;
                }

                if (!schema.hasHeader && !headerWritten) {
                    String[] generated = new String[fields.length];
                    for (int i = 0; i < fields.length; i++) generated[i] = "col" + (i + 1);
                    bw.write(InputTable.join(append(generated, ANNOTATION_HEADER), schema.delimiter));
                    bw.newLine();
                    headerWritten = true;
                }

                rows++;
                try {
                    if (schema.chromosomeIndex >= fields.length ||
                        schema.startIndex >= fields.length ||
                        schema.endIndex >= fields.length) {
                        throw new IllegalArgumentException("row has fewer columns than the detected schema");
                    }

                    String chr = InputTable.stripQuotes(fields[schema.chromosomeIndex]);
                    int rawStart = Integer.parseInt(InputTable.stripQuotes(fields[schema.startIndex]));
                    int rawEnd = Integer.parseInt(InputTable.stripQuotes(fields[schema.endIndex]));

                    int queryStart, queryEnd;
                    if (schema.coordinates == InputTable.CoordinateMode.BED) {
                        queryStart = rawStart + 1;
                        queryEnd = rawEnd;
                    } else {
                        queryStart = rawStart;
                        queryEnd = rawEnd;
                    }
                    if (queryEnd < queryStart) {
                        int t = queryStart;
                        queryStart = queryEnd;
                        queryEnd = t;
                    }

                    String strand = null;
                    if (schema.strandIndex >= 0 && schema.strandIndex < fields.length) {
                        strand = InputTable.normalizeStrandToken(
                                InputTable.stripQuotes(fields[schema.strandIndex]));
                    }

                    RnaAnnotationEngine.Result result = engine.annotate(chr, queryStart, queryEnd, strand);
                    String[] annotation = new String[]{
                            result.geneId, result.geneName, result.transcriptIds,
                            result.geneType, result.transcriptType,
                            result.geneStrand, result.location
                    };

                    bw.write(InputTable.join(append(fields, annotation), schema.delimiter));
                    bw.newLine();
                    written++;
                } catch (Exception rowError) {
                    invalid++;
                    if (!cli.quiet && invalid <= 5) {
                        System.err.println("WARNING: skipped malformed row " + rows + ": " + rowError.getMessage());
                    }
                }
            }
        }

        if (!cli.quiet) {
            System.err.println("Done. Annotated rows: " + written + "; skipped rows: " + invalid);
            System.err.println("Output: " + cli.output);
        }
    }

    private static Cli parse(String[] args) {
        if (args.length == 0) {
            usage();
            return null;
        }

        for (String arg : args) {
            if ("--help".equals(arg) || "-h".equals(arg)) {
                usage();
                return null;
            }
            if ("--version".equals(arg)) {
                System.out.println("PeakLocus " + VERSION);
                return null;
            }
        }

        Cli cli = new Cli();
        for (int i = 0; i < args.length; i++) {
            String arg = args[i];
            if ("--input".equals(arg) || "-i".equals(arg)) cli.input = Paths.get(require(args, ++i, arg));
            else if ("--gtf".equals(arg) || "-g".equals(arg)) cli.gtf = Paths.get(require(args, ++i, arg));
            else if ("--output".equals(arg) || "-o".equals(arg)) cli.output = Paths.get(require(args, ++i, arg));
            else if ("--strand-column".equals(arg)) cli.strandColumn = require(args, ++i, arg);
            else if ("--chrom-column".equals(arg)) cli.chromosomeColumn = require(args, ++i, arg);
            else if ("--start-column".equals(arg)) cli.startColumn = require(args, ++i, arg);
            else if ("--end-column".equals(arg)) cli.endColumn = require(args, ++i, arg);
            else if ("--coordinates".equals(arg)) {
                String value = require(args, ++i, arg).toLowerCase(Locale.ROOT);
                if ("auto".equals(value)) cli.coordinates = InputTable.CoordinateMode.AUTO;
                else if ("bed".equals(value) || "0-based".equals(value)) cli.coordinates = InputTable.CoordinateMode.BED;
                else if ("one-based".equals(value) || "1-based".equals(value)) cli.coordinates = InputTable.CoordinateMode.ONE_BASED;
                else throw new IllegalArgumentException("Unknown --coordinates: " + value);
            }
            else if ("--quiet".equals(arg)) cli.quiet = true;
            else throw new IllegalArgumentException("Unknown argument: " + arg);
        }

        if (cli.input == null || cli.gtf == null) {
            throw new IllegalArgumentException("--input and --gtf are required. Use --help for usage.");
        }
        if (!Files.isRegularFile(cli.input)) throw new IllegalArgumentException("Input not found: " + cli.input);
        if (!Files.isRegularFile(cli.gtf)) throw new IllegalArgumentException("GTF not found: " + cli.gtf);
        return cli;
    }

    private static String require(String[] args, int index, String flag) {
        if (index >= args.length) throw new IllegalArgumentException("Missing value after " + flag);
        return args[index];
    }

    private static String[] append(String[] left, String[] right) {
        String[] out = Arrays.copyOf(left, left.length + right.length);
        System.arraycopy(right, 0, out, left.length, right.length);
        return out;
    }

    private static Path defaultOutput(Path input, InputTable.Delimiter delimiter) {
        String name = input.getFileName().toString();
        if (name.endsWith(".gz")) name = name.substring(0, name.length() - 3);
        int dot = name.lastIndexOf('.');
        String base = dot > 0 ? name.substring(0, dot) : name;
        String extension = delimiter == InputTable.Delimiter.COMMA ? ".peaklocus.csv" : ".peaklocus.tsv";
        Path parent = input.getParent();
        return parent == null ? Paths.get(base + extension) : parent.resolve(base + extension);
    }

    private static void usage() {
        System.out.println(
            "PeakLocus v" + VERSION + "\n" +
            "RNA-centric genomic peak and interval annotation.\n\n" +
            "PeakLocus maps RNA-associated genomic intervals to genes, transcripts,\n" +
            "gene/transcript biotypes, strand, and RNA-region context:\n" +
            "3'UTR, 5'UTR, CDS, exon, intron, or N/A.\n\n" +
            "Usage:\n" +
            "  java -jar peaklocus.jar -i INPUT -g ANNOTATION.gtf [-o OUTPUT] [options]\n\n" +
            "Required:\n" +
            "  -i, --input FILE                  RNA peak/interval table\n" +
            "  -g, --gtf FILE                    GTF annotation (plain or .gz)\n\n" +
            "Automatic input detection:\n" +
            "  --strand-column auto|none|N|NAME  Default: auto\n" +
            "  --chrom-column auto|N|NAME        Default: auto\n" +
            "  --start-column auto|N|NAME        Default: auto\n" +
            "  --end-column auto|N|NAME          Default: auto\n" +
            "  --coordinates auto|bed|one-based  Default: auto\n\n" +
            "Other:\n" +
            "  -o, --output FILE                 Default: INPUT.peaklocus.tsv\n" +
            "  --quiet                           Suppress progress messages\n" +
            "  --version                         Print version\n" +
            "  -h, --help                        Show this help\n\n" +
            "Column numbers supplied on the command line are 1-based.\n");
    }
}
