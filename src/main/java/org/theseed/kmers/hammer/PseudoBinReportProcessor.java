/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;

import org.slf4j.LoggerFactory;
import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.theseed.counters.WeightMap;
import org.theseed.io.TabbedLineReader;
import org.theseed.utils.BaseReportProcessor;
import org.theseed.utils.ParseFailureException;


/**
 * This command produces a simulated bin report from binning results.  The input directory is a master sample
 * group that has already been binned.  That is, it is a directory whose subdirectories are samples that have
 * been binned in place.  Each sample subdirectory should contain the "report.tbl" file produced at the end of
 * binning.  This report is used to generate a simulated bin report for comparison with the hammer bin reports.
 *
 * The positional parameters should be the name of the input sample master directory and the name of the
 * appropriate repgen list file.  The list file is used to convert the genome IDs in the "report.tbl" file
 * to the representative genomes used in the hammer bin reports.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output file for bin report (if not STDOUT)
 *
 * @author Bruce Parrello
 *
 */
public class PseudoBinReportProcessor extends BaseReportProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(PseudoBinReportProcessor.class);
    /** map of genome IDs to representative genome IDs */
    private Map<String, String> repGenMap;
    /** map of representative genome IDs to names */
    private Map<String, String> repNameMap;
    /** list of sample directories */
    private File[] sampleDirs;
    /** file filter for sample directories */
    private FileFilter SAMPLE_DIR_FILTER = new FileFilter() {
        @Override
        public boolean accept(File pathname) {
            boolean retVal = pathname.isDirectory();
            if (retVal) {
                File reportFile = new File(pathname, "report.tbl");
                retVal = reportFile.isFile();
            }
            return retVal;
        }
    };

    // COMMAND-LINE OPTIONS

    @Argument(index = 0, metaVar = "inDir", usage = "input sample master directory", required = true)
    private File inDir;

    @Argument(index = 1, metaVar = "repgen.list.tbl", usage = "repgen list file used to create hammers")
    private File repListFile;

    @Override
    protected void setReporterDefaults() {
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        // Verify that the input directory exists.
        if (! this.inDir.isDirectory())
            throw new FileNotFoundException("Input master directory " + this.inDir + " is not found or unreadable.");
        // Insure that it has binned samples in it.
        this.sampleDirs = this.inDir.listFiles(SAMPLE_DIR_FILTER);
        if (this.sampleDirs.length == 0)
            throw new IOException("No binned samples found in " + this.inDir + ".");
        log.info("{} binned samples found in {}.", this.sampleDirs.length, this.inDir);
        // Verify and process the repgen list file.
        if (! this.repListFile.canRead())
            throw new FileNotFoundException("Repgen list file " + this.repListFile + " is not found or unreadable.");
        this.readRepListFile();
    }

    /**
     * Read in the repgen list file.  We need a map of genome IDs to representative genome IDs, and a map
     * of representative genome IDs to genome names.  Since the names are only provided for the actual genome
     * IDs, we storea genome name if it is its own representative.
     *
     * @throws IOException
     */
    private void readRepListFile() throws IOException {
        // Create the maps.
        this.repGenMap = new HashMap<String, String>(1000);
        this.repNameMap = new HashMap<String, String>(1000);
        // Process the file to fill the maps.
        try (TabbedLineReader repListStream = new TabbedLineReader(this.repListFile)) {
            int gCol = repListStream.findField("genome_id");
            int nameCol = repListStream.findField("genome_name");
            int repCol = repListStream.findField("rep_id");
            log.info("Loading genome reprresentation data from {}.", this.repListFile);
            for (var line : repListStream) {
                String genomeId = line.get(gCol);
                String repId = line.get(repCol);
                this.repGenMap.put(genomeId, repId);
                if (genomeId.contentEquals(repId))
                    this.repNameMap.put(genomeId, line.get(nameCol));
            }
        }
        log.info("{} genomes loaded, {} representatives found.", this.repGenMap.size(), this.repNameMap.size());
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // Compute the optimal hash size.
        int hashSize = this.repNameMap.size() * 4 / 3 + 1;
        // Write the report header.
        writer.println("sample_id\trepgen_id\trepgen_name\tcount");
        // Now we loop through the samples, using the bin report to compute weights.
        for (File sampleDir : this.sampleDirs) {
            String sampleId = sampleDir.getName();
            log.info("Processing sample {}.", sampleId);
            WeightMap weights = new WeightMap(hashSize);
            // Read the report file.
            File reportFile = new File(sampleDir, "report.tbl");
            try (TabbedLineReader reportStream = new TabbedLineReader(reportFile)) {
                int refCol = reportStream.findField("ref_genome_id");
                int covgCol = reportStream.findField("coverage");
                int sizeCol = reportStream.findField("dna_size");
                for (var line : reportStream) {
                    // Compute the scale factor for this bin.
                    double scale = line.getInt(sizeCol) * line.getDouble(covgCol) / 170;
                    // Break out the reference genome IDs.
                    String[] refGenomes = StringUtils.split(line.get(refCol), ", ");
                    for (String refGenome : refGenomes) {
                        // Get the representative for this genome.
                        String repId = this.repGenMap.get(refGenome);
                        if (repId == null)
                            log.warn("ERROR: genome {} is not in the repgen map.", refGenome);
                        else {
                            // Record it in the scoring map.
                            weights.count(repId, scale);
                        }
                    }
                }
                log.info("Writing scores for {}.", sampleId);
                var scores = weights.sortedCounts();
                for (var score : scores) {
                    String repId = score.getKey();
                    double count = score.getCount();
                    String repName = this.repNameMap.getOrDefault(repId, "<unknown>");
                    writer.println(sampleId + "\t" + repId + "\t" + repName + "\t" + count);
                }
            }
        }
    }

}
