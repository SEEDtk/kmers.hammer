/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Map;

import org.apache.commons.io.filefilter.SuffixFileFilter;
import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.counters.WeightMap;
import org.theseed.io.TabbedLineReader;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.Sequence;
import org.theseed.utils.BaseReportProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command reads a synthetic-sample FASTA file and converts it into a bin report.  It is used to compare
 * with bin reports produced by other methods (hammers, binning, Kraken).
 *
 * The positional parameters are the name of the input directory containing the synthetic samples and the name of the
 * repXXX.stats.tbl file containing the repgen genome IDs and names.  The samples should all have the filename suffix
 * ".fa" and the base part of the file name should be the sample name.  Each contig should have a tab-delimited comment
 * with a repgen ID in the second field.  The score for each repgen will be the sum of the relevant contig lengths
 * divided by a constant (currently 100K).
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output file for report (if not STDOUT)
 *
 * @author Bruce Parrello
 *
 */
public class SynthBinReportProcessor extends BaseReportProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(SynthBinReportProcessor.class);
    /** array of FASTA files to read */
    private File[] inFiles;
    /** map of repgen IDs to names */
    private Map<String, String> repNameMap;
    /** scale factor for converting contig lengths to scores */
    private final static double SCALE_FACTOR = 100000.0;

    // COMMAND-LINE OPTIONS

    /** input directory */
    @Argument(index = 0, metaVar = "inDir", usage = "input directory containing sample FASTA files", required = true)
    private File inDir;

    /** repgen stats file */
    @Argument(index = 1, metaVar = "repXXX.stats.tbl", usage = "repgen stats file containing representative genome IDs and names",
            required = true)
    private File repStatsFile;

    @Override
    protected void setReporterDefaults() {
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        // Read in the repgen stats file to create the name map.
        if (! this.repStatsFile.canRead())
            throw new FileNotFoundException("Repgen stats file " + this.repStatsFile + " is not found or unreadable.");
        this.repNameMap = TabbedLineReader.readMap(this.repStatsFile, "rep_id", "rep_name");
        log.info("{} repgen genomes found.", this.repNameMap.size());
        // Get the list of files in the input directory.
        if (! this.inDir.isDirectory())
            throw new FileNotFoundException("Input directory " + this.inDir + " is not found or unreadable.");
        this.inFiles = this.inDir.listFiles((FileFilter) new SuffixFileFilter(".fa"));
        log.info("{} sample FASTA files found in {}.", this.inFiles.length, this.inDir);
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // Write the output header.
        writer.println("sample_id\trepgen_id\trepgen_name\tcount");
        // This map will store the scores.
        WeightMap scoreMap = new WeightMap(this.repNameMap.size() * 4 / 3 + 1);
        // Loop through the FASTA files.
        for (File inFile : this.inFiles) {
            String sampleId = StringUtils.substringBeforeLast(inFile.getName(), ".");
            log.info("Processing sample {} in file {}.", sampleId, inFile);
            // Clear the score map.
            scoreMap.deleteAll();
            // Read the contigs in the file.
            int count = 0;
            int errors = 0;
            try (var contigStream = new FastaInputStream(inFile)) {
                for (Sequence contig : contigStream) {
                    // Get the contig length and convert it to a score.
                    double score = contig.length() / SCALE_FACTOR;
                    String[] commentParts = StringUtils.split(contig.getComment(), '\t');
                    String rep_id = commentParts[1];
                    String rep_name = this.repNameMap.get(rep_id);
                    count++;
                    if (rep_name == null) {
                        log.error("Invalid repgen ID {} in contig {} of sample {}.", rep_id, contig.getLabel(), sampleId);
                        errors++;
                    } else
                        scoreMap.count(rep_id, score);
                    if (count % 1000 == 0)
                        log.info("{} contigs processed in sample {}. {} errors.", count, sampleId, errors);
                }
                log.info("{} contigs processed in sample {}. {} errors.", count, sampleId, errors);
                // Output the weight map for this sample.
                for (var repEntry : scoreMap.counts()) {
                    String rep_id = repEntry.getKey();
                    String rep_name = this.repNameMap.get(rep_id);
                    double score = repEntry.getCount();
                    writer.println(sampleId + "\t" + rep_id + "\t" + rep_name + "\t" + Double.toString(score));
                }
                writer.flush();
            }
        }
    }

}
