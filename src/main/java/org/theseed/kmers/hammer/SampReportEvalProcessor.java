/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UncheckedIOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.io.TabbedLineReader;
import org.theseed.utils.BaseMultiReportProcessor;
import org.theseed.utils.StringPair;
import org.theseed.reports.eval.SampReportEvalReporter;
import org.theseed.reports.eval.SampReportEvalReporter.SampleDescriptor;
import org.theseed.stats.WeightMap;

/**
 *
 * This command analyzes the results of one or more hammer tests against FASTQ samples and compares the results.
 * The bin-report output files are specified as positional parameters on the command line.  The standard input
 * should contain a report characterizing the samples, all of which are expected to be single-organism read sets.
 * For each sample, the report should have the sample ID in a column named "sample_id", the organism ID and name in
 * "genome_id" and "genome_name", and the expected representative ID in a column named by the command-line options.
 * The input file can contain extra samples:  these will be ignored. Multiple output reports will be produced in
 * a single output directory.
 *
 * Some output reports can make use of genome distance information.  In this case, a distance file must be provided.  The
 * genome IDs must be in columns named "id1" and "id2".  The distance value itself is specified by a parameter,
 * but defaults to the last column.
 *
 * The input reports will be identified by base name with the extension removed.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -D	output directory name (default SampEval in the current directory)
 * -i	input file for sample characterization (if not STDIN)
 * -c	index (1-based) or name of input column name for expected representative (default "rep_id")
 *
 * --distFile	name of a file containing genome distances, from the genome.distance methods command (optional)
 * --distCol	index (1-based) or name of column containing the distance to use in the distance file (default "0")
 * --roles		tab-delimited file with headers containing a list of role IDs in the first column, for roles of special interest
 * --clear		erase the output directory before processing
 *
 * @author Bruce Parrello
 *
 */
public class SampReportEvalProcessor extends BaseMultiReportProcessor implements SampReportEvalReporter.IParms {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(SampReportEvalProcessor.class);
    /** input column index for sample ID */
    private int sampleIdIdx;
    /** input column index for genome ID */
    private int genomeIdIdx;
    /** input column index for genome name */
    private int genomeNameIdx;
    /** input column index for expected representative */
    private int repIdIdx;
    /** sample description map */
    private Map<String, Descriptor> sampMap;
    /** output reporters */
    private List<Report> reports;
    /** distance table */
    private Map<StringPair, Double> distTable;
    /** role set */
    private Set<String> roleSet;

    // COMMAND-LINE OPTIONS

    /** expected-representative input column name */
    @Option(name = "--col", aliases = { "-c" }, metaVar = "rep100", usage = "index (1-based) or name of the input column for the expected representative")
    private String repCol;

    /** input file (if not STDIN) containing sample data */
    @Option(name = "--input", aliases = { "-i" }, metaVar = "inFile.tbl", usage = "input file (if not STDIN)")
    private File inFile;

    /** distance file to use for distance-based reports */
    @Option(name = "--distFile", aliases = { "--dists" }, metaVar = "genomes.dists.tbl", usage = "name of an optional file containing genome distances")
    private File distFile;

    /** column containing the distances in the distance file */
    @Option(name = "--distCol", metaVar = "ANI", usage = "index (1-based) or name of distance column in distance file")
    private String distCol;

    /** file containing roles of interest */
    @Option(name = "--roles", metaVar = "roles.tbl", usage = "optional tab-delimited file of role IDs for roles of interest")
    private File roleFile;

   /** list of input sample bin reports */
    @Argument(index = 0, metaVar = "reportFile1 reportFile2 ...", usage = "names of files containing sample bin reports to evaluate (wildcards allowed)")
    private List<File> reportFiles;

    /**
     * This class describes the data we need for a sample from the input file.
     */
    public class Descriptor extends SampReportEvalReporter.SampleDescriptor {

        /**
         * Create the sample descriptor from an input line.
         *
         * @param line	input line for the sample
         */
        protected Descriptor(TabbedLineReader.Line line) {
            this.setSampleId(line.get(SampReportEvalProcessor.this.sampleIdIdx));
            this.setGenomeId(line.get(SampReportEvalProcessor.this.genomeIdIdx));
            this.setGenomeName(line.get(SampReportEvalProcessor.this.genomeNameIdx));
            this.setRepId(line.get(SampReportEvalProcessor.this.repIdIdx));
        }

    }

    /**
     * This class contains a reporter object and its associated report writer.
     */
    protected class Report {

        /** reporter object */
        private SampReportEvalReporter reporter;
        /** output print writer */
        private PrintWriter writer;

        /**
         * Construct the reporter for a specified report type
         *
         * @param type	type of report
         *
         * @throws IOException
         */
        protected Report(SampReportEvalReporter.Type type) throws IOException {
            this.writer = SampReportEvalProcessor.this.openReport(StringUtils.lowerCase(type.name()) + ".tbl");
            this.reporter = type.create(SampReportEvalProcessor.this, this.writer);
        }

        /**
         * Begin processing a new bin report.
         *
         * @param reportName	name of bin report
         */
        public void openFile(String name) {
            this.reporter.openFile(name);
        }

        /**
         * Record hits against a sample.
         *
         * @param desc			relevant sample descriptor
         * @param repId			ID of the repgen hit
         * @param repName		name of the repgen hit
         * @param count			weighted hit count
         * @param roleCount 	number of roles hit
         * @param roleCounts 	map of number of hits per each role
         *
         */
        public void recordHits(SampleDescriptor desc, String repId, String repName, double count, int roleCount,
                WeightMap roleCounts) {
            try {
                this.reporter.recordHits(desc, repId, repName, count, roleCount, roleCounts);
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }

        }

        /**
         * Finish processing a bin report.
         */
        public void closeFile() {
            try {
                this.reporter.closeFile();
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        }

        /**
         * Finish the entire report.
         */
        public void closeReport() {
            this.reporter.closeReport();
        }

        /**
         * Close the report writer.
         */
        public void close() {
            this.writer.close();
        }

    }

    @Override
    protected File setDefaultOutputDir(File curDir) {
        return new File(System.getProperty("user.dir"), "SampEval");
    }

    @Override
    protected void setMultiReportDefaults() {
        this.repCol = "rep_id";
        this.reportFiles = new ArrayList<File>();
        this.distCol = "0";
        this.roleFile = null;
    }

    @Override
    protected void validateMultiReportParms() throws IOException, ParseFailureException {
        // Here we validate the input report files.
        log.info("Checking input reports.");
        if (this.reportFiles.isEmpty())
            throw new ParseFailureException("No input report files specified.");
        for (File reportFile : this.reportFiles) {
            if (! reportFile.canRead())
                throw new FileNotFoundException("Input file " + reportFile + " is not found or unreadable.");
            else
                log.debug("Input file recorded: {}.", reportFile);
        }
        log.info("{} input report files specified.", this.reportFiles.size());
        // Read in the sample characterizations.
        try (var inputStream = TabbedLineReader.openInput(this.inFile)) {
            // Find all the input columns.
            this.genomeIdIdx = inputStream.findField("genome_id");
            this.genomeNameIdx = inputStream.findField("genome_name");
            this.sampleIdIdx = inputStream.findField("sample_id");
            this.repIdIdx = inputStream.findField(this.repCol);
            // Create and fill the sample map.
            this.sampMap = new HashMap<String, Descriptor>();
            for (var line : inputStream) {
                Descriptor desc = this.new Descriptor(line);
                this.sampMap.put(desc.getSampleId(), desc);
            }
            log.info("{} samples found in input.", this.sampMap.size());
        }
        // Check for a distance file.
        this.distTable = new HashMap<StringPair, Double>();
        if (this.distFile != null) {
            log.info("Loading distances from {}.", this.distFile);
            try (TabbedLineReader distStream = new TabbedLineReader(this.distFile)) {
                int id1ColIdx = distStream.findField("id1");
                int id2ColIdx = distStream.findField("id2");
                int distColIdx = distStream.findField(this.distCol);
                for (var line : distStream) {
                    String g1 = line.get(id1ColIdx);
                    String g2 = line.get(id2ColIdx);
                    double dist = line.getDouble(distColIdx);
                    StringPair gPair = new StringPair(g1, g2);
                    this.distTable.put(gPair, dist);
                }
                log.info("{} distances stored.", this.distTable.size());
            }
        }
        // Check for a role file.
        if (this.roleFile == null)
            this.roleSet = Collections.emptySet();
        else {
            // Note we need to sort the role IDs when we read them in.
            this.roleSet = new TreeSet<String>(TabbedLineReader.readSet(this.roleFile, "1"));
            log.info("{} roles of interest loaded from {}.", this.roleSet.size(), this.roleFile);
        }
    }

    @Override
    protected void runMultiReports() throws Exception {
        try {
            // Here we create our reports and the associated output files.
            this.reports = new ArrayList<Report>(SampReportEvalReporter.Type.values().length);
            for (SampReportEvalReporter.Type type : SampReportEvalReporter.Type.values()) {
                Report newReport = this.new Report(type);
                this.reports.add(newReport);
            }
            log.info("{} reports initialized in {}.", this.reports.size(), this.getOutDir());
            // Now we loop through the report files.
            for (File reportFile : this.reportFiles) {
                String reportName = StringUtils.substringBeforeLast(reportFile.getName(), ".");
                log.info("Processing report file {} with name {}.", reportFile, reportName);
                this.reports.stream().forEach(x -> x.openFile(reportName));
                try (TabbedLineReader reportStream = new TabbedLineReader(reportFile)) {
                    // Create a count map for the role columns.
                    String[] labels = reportStream.getLabels();
                    WeightMap roleCounts = new WeightMap();
                    // Loop through the report file.
                    for (var line : reportStream) {
                        String sampleId = line.get(0);
                        Descriptor desc = this.sampMap.get(sampleId);
                        if (desc == null)
                            throw new IOException("Sample " + sampleId + " is not in the sample input file.");
                        String repId = line.get(1);
                        String repName = line.get(2);
                        double count = line.getDouble(3);
                        int roleCount = line.getInt(4);
                        for (int i = 5; i < reportStream.size(); i++)
                            roleCounts.setCount(labels[i], line.getDouble(i));
                        this.reports.stream().forEach(x -> x.recordHits(desc, repId, repName, count, roleCount, roleCounts));
                    }
                    // All done, finish the reporting for this sample.
                    this.reports.stream().forEach(x -> x.closeFile());
                }
            }
            // Finish the full report.
            this.reports.stream().forEach(x -> x.closeReport());
        } finally {
            // Close the output streams.
            this.reports.stream().forEach(x -> x.close());
        }
    }


    @Override
    public SampReportEvalReporter.SampleDescriptor getDescriptor(String sampleId) {
        return this.sampMap.get(sampleId);
    }

    @Override
    public double getDistance(String genome1, String genome2) {
        StringPair gPair = new StringPair(genome1, genome2);
        return this.distTable.getOrDefault(gPair, Double.NaN);
    }

    @Override
    public Set<String> getRoleSet() {
        return this.roleSet;
    }


}
