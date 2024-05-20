/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.io.TabbedLineReader;
import org.theseed.reports.eval.SampReportEvalReporter;
import org.theseed.reports.eval.SampReportEvalReporter.SampleDescriptor;
import org.theseed.utils.BasePipeProcessor;

/**
 * This command analyzes the resullts of one or more hammer tests against FASTQ samples and compares the results.
 * The bin-report output files are specified as positional parameters on the command line.  The standard input
 * should contain a report characterizing the samples, all of which are expected to be single-organism read sets.
 * For each sample, the report should have the sample ID in a column named "sample_id", the organism ID and name
 * in "genome_id" and "genome_name", and the expected representative ID in a column named by the command-line options.
 * The input file can contain extra samples:  these will be ignored.
 *
 * The reports will be identified by base name with the extension removed.  Multiple report types are supported.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -i	input file for sample characterization (if not STDIN)
 * -o	output file for report (if not STDOUT)
 * -c	index (1-based) or name of input column name for expected representative (default "rep_id")
 *
 * --format		type of report to output (default QUALITY)
 *
 * @author Bruce Parrello
 *
 */
public class SampReportEvalProcessor extends BasePipeProcessor implements SampReportEvalReporter.IParms {

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
    /** output reporter */
    private SampReportEvalReporter reporter;

    // COMMAND-LINE OPTIONS

    /** expected-representative input column name */
    @Option(name = "--col", aliases = { "-c" }, metaVar = "rep100", usage = "index (1-based) or name of the input column for the expected representative")
    private String repCol;

    /** output format for the report */
    @Option(name = "--format", metaVar = "output format for the report")
    private SampReportEvalReporter.Type reportType;

    /** list of input sample bin reports */
    @Argument(index = 0, metaVar = "reportFile1 reportFile2 ...", usage = "files containing sample bin reports to evaluate")
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


    @Override
    protected void setPipeDefaults() {
        this.repCol = "rep_id";
        this.reportFiles = new ArrayList<File>();
        this.reportType = SampReportEvalReporter.Type.QUALITY;
    }

    @Override
    protected void validatePipeParms() throws IOException, ParseFailureException {
        // Here we validate the input report files.
        if (this.reportFiles.isEmpty())
            throw new ParseFailureException("No input report files specified.");
        for (File reportFile : this.reportFiles) {
            if (! reportFile.canRead())
                throw new FileNotFoundException("Input file " + reportFile + " is not found or unreadable.");
        }
        log.info("{} input report files specified.", this.reportFiles.size());
    }

    @Override
    protected void validatePipeInput(TabbedLineReader inputStream) throws IOException {
        // Find all the input columns.
        this.genomeIdIdx = inputStream.findField("genome_id");
        this.genomeNameIdx = inputStream.findField("genome_name");
        this.sampleIdIdx = inputStream.findField("sample_id");
        this.repIdIdx = inputStream.findField(this.repCol);
    }

    @Override
    protected void runPipeline(TabbedLineReader inputStream, PrintWriter writer) throws Exception {
        // Create the output reporter.
        this.reporter = this.reportType.create(this, writer);
        // Read in the sample characterizations.
        log.info("Reading sample characterizations.");
        this.sampMap = new HashMap<String, Descriptor>();
        for (var line : inputStream) {
            Descriptor desc = this.new Descriptor(line);
            this.sampMap.put(desc.getSampleId(), desc);
        }
        log.info("{} samples found in input file.", this.sampMap.size());
        // Now we loop through the report files.
        for (File reportFile : this.reportFiles) {
            String reportName = StringUtils.substringBeforeLast(reportFile.getName(), ".");
            log.info("Processing report file {} with name {}.", reportFile, reportName);
            this.reporter.openFile(reportName);
            try (TabbedLineReader reportStream = new TabbedLineReader(reportFile)) {
                // Loop through the report file.
                for (var line : reportStream) {
                    String sampleId = line.get(0);
                    Descriptor desc = this.sampMap.get(sampleId);
                    String repId = line.get(1);
                    String repName = line.get(2);
                    double count = line.getDouble(3);
                    this.reporter.recordHits(desc, repId, repName, count);
                }
                // All done, finish the reporting for this sample.
                this.reporter.closeFile();
            }
        }
        // Finish the full report.
        this.reporter.closeReport();
    }

    @Override
    public SampleDescriptor getDescriptor(String sampleId) {
        return this.sampMap.get(sampleId);
    }

}
