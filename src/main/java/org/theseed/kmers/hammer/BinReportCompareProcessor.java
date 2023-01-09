/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UncheckedIOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.apache.commons.math3.util.ResizableDoubleArray;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.counters.WeightMap;
import org.theseed.io.TabbedLineReader;
import org.theseed.utils.BaseReportProcessor;
import org.theseed.utils.IntegerList;
import org.theseed.utils.ParseFailureException;

/**
 * This command compares bin reports.  We will load each report into a weight map and compare the weight maps on a sample-by
 * sample basis.  Given a sample and two reports, we first do spearman and pearson correlations.  Then, we do a top-N comparison
 * for various values of N.  The comparison will only be performed if both maps have at least N groupings present, and will
 * essentially be a similarity ratio-- number of groupings in common over size of set.
 *
 * The positional parameters are the names of the bin report files to compare.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output file for report (if not STDOUT)
 *
 * --N		comma-delimited list of N-values to use for top-N comparisons (default 5,10,15,20)
 *
 * @author Bruce Parrello
 *
 */
public class BinReportCompareProcessor extends BaseReportProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(BinReportCompareProcessor.class);
    /** engine for spearman's correlation */
    private SpearmansCorrelation spearmanEngine;
    /** engine for pearson's coefficient */
    private PearsonsCorrelation pearsonEngine;
    /** list of N-values */
    private IntegerList nValueList;
    /** list of (sampleId -> weightMap) maps for each file */
    private List<Map<String, WeightMap>> fileContents;

    // COMMAND-LINE OPTIONS

    /** list of N-values to use for top-N stuff */
    @Option(name = "--N", metaVar = "4,8,12,16,20", usage = "comma-delimited list of N-values to use for top-N similarities")
    private String nValueString;

    /** file names of bin reports to compare */
    @Argument(index = 0, metaVar = "file1.tbl file2.tbl ...", usage = "file name(s) of bin report(s) to compare", required = true)
    private List<File> reportFiles;

    @Override
    protected void setReporterDefaults() {
        this.nValueString = "5,10,15,20";
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        // Validate the n-value list.
        try {
            this.nValueList = new IntegerList(this.nValueString);
        } catch (NumberFormatException e) {
            throw new ParseFailureException("Invalid N-value list \"" + this.nValueString + "\": must be valid integers separated by commas.");
        }
        for (int n : this.nValueList) {
            if (n <= 0)
                throw new ParseFailureException("Cannot use " + Integer.toString(n) + " in the N-value list because it is less than 1.");
        }
        if (this.nValueList.size() < 1)
            throw new ParseFailureException("Must be at least one value in N-value list.");
        // Now check all the bin reports.
        for (File reportFile : this.reportFiles) {
            if (! reportFile.canRead())
                throw new FileNotFoundException("Input file " + reportFile + " is not found or unreadable.");
        }
        if (this.reportFiles.size() < 2)
            throw new ParseFailureException("Must be at least 2 input files specified.");
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // Create the correlation engines.
        this.spearmanEngine = new SpearmansCorrelation();
        this.pearsonEngine = new PearsonsCorrelation();
        // Write the header line.
        writer.print("file1\tfile2\tsample_id\tspearman\tpearson\tall_sim");
        for (int n : this.nValueList)
            writer.format("\ttop%d_sim", n);
        writer.println();
        // Create the list of file-content maps.
        final int nFiles = this.reportFiles.size();
        this.fileContents = this.reportFiles.stream().map(x -> this.readReport(x)).collect(Collectors.toList());
        // Insure each file has the same list of samples.
        Set<String> samples = this.fileContents.get(0).keySet();
        for (var fileContent : this.fileContents.subList(1, nFiles)) {
            if (! samples.equals(fileContent.keySet()))
                throw new IOException("All report files must contain the same set of samples.");
        }
        // Process each sample in each pair of files.
        for (int i1 = 0; i1 < nFiles; i1++) {
            for (int i2 = i1 + 1; i2 < nFiles; i2++) {
                // Get the content and name for each file.
                var file1Map = this.fileContents.get(i1);
                var file2Map = this.fileContents.get(i2);
                String file1Name = this.reportFiles.get(i1).getName();
                String file2Name = this.reportFiles.get(i2).getName();
                // Loop through the samples.
                for (String sampleId : samples) {
                    log.info("Comparing files {} and {} for {}.", file1Name, file2Name, sampleId);
                    // Get the weight maps for this sample.
                    WeightMap scores1 = file1Map.get(sampleId);
                    WeightMap scores2 = file2Map.get(sampleId);
                    // Ask for the metrics on these weight maps.
                    double[] metrics = this.compareMaps(scores1, scores2);
                    // Write out the results.
                    writer.println(file1Name + "\t" + file2Name + "\t" + sampleId + "\t"
                            + Arrays.stream(metrics).mapToObj(x -> Double.toString(x)).collect(Collectors.joining("\t")));
                }
            }
        }
    }

    /**
     * Read in a report file and return a map of weightmaps for each sample.
     *
     * @param reportFile	bin-report file to read
     *
     * @return a map from each sample ID in the bin report to its weight map based on the counts in the file
     */
    private Map<String, WeightMap> readReport(File reportFile) {
        // Create the return map.
        Map<String, WeightMap> retVal = new HashMap<String, WeightMap>();
        // This will count the score lines read.
        int scoreCount = 0;
        // Set up the file.
        try (TabbedLineReader reportStream = new TabbedLineReader(reportFile)) {
            int sampleIdCol = reportStream.findField("sample_id");
            int repIdCol = reportStream.findField("repgen_id");
            int scoreCol = reportStream.findField("count");
            // Loop through the file, building the weight maps.
            for (var line : reportStream) {
                String sampleId = line.get(sampleIdCol);
                String repId = line.get(repIdCol);
                double score = line.getDouble(scoreCol);
                WeightMap scores = retVal.computeIfAbsent(sampleId, x -> new WeightMap());
                scores.count(repId, score);
                scoreCount++;
            }
        } catch (IOException e) {
            // Convert an IO error to an unchecked exception so we can use this method in a stream.
            throw new UncheckedIOException(e);
        }
        log.info("{} samples and {} scores read from {}.", retVal.size(), scoreCount, reportFile);
        return retVal;
    }

    /**
     * Compute the array of metrics for a pair of score maps.
     *
     * @param scores1	first score map
     * @param scores2	second score map
     *
     * @return an array of comparison measures, starting with the three all-map measures and then the top-N measures
     */
    private double[] compareMaps(WeightMap scores1, WeightMap scores2) {
        // First we need to compute the correlations.  We use key-sorted counts for this.
        var counts1 = scores1.keyedCounts();
        var counts2 = scores2.keyedCounts();
        // Do a classic merge to put the counts in a pair of arrays.
        int i1 = 0;
        var array1 = new ResizableDoubleArray(counts1.size());
        int i2 = 0;
        var array2 = new ResizableDoubleArray(counts2.size());
        // We stop at the end of the first array, so we can only fall off the end of the second inside the loop.
        while (i1 < counts1.size()) {
            // Compare the left and right keys.  If there is no right key, the left key is less (negative).
            int comparison = (i2 >= counts2.size() ? -1 : counts1.get(i1).getKey().compareTo(counts2.get(i2).getKey()));
            if (comparison <= 0) {
                array1.addElement(counts1.get(i1).getCount());
                i1++;
            } else
                array1.addElement(0.0);
            if (comparison >= 0) {
                array2.addElement(counts2.get(i2).getCount());
                i2++;
            } else
                array2.addElement(0.0);
        }
        // Here there might be leftovers on the right.
        while (i2 < counts2.size()) {
            array1.addElement(0.0);
            array2.addElement(counts2.get(i2).getCount());
            i2++;
        }
        double[] double1 = array1.getElements();
        double[] double2 = array2.getElements();
        // Now we have two parallel arrays, each containing a position for every single representative in either map.
        // The size of either array is the denominator for the all-similarity metric.  Allocate an array for the
        // output metrics.
        double[] retVal = new double[3 + this.nValueList.size()];
        retVal[0] = this.spearmanEngine.correlation(double1, double2);
        retVal[1] = this.pearsonEngine.correlation(double1, double2);
        // Compute the all-similarity ratio.
        retVal[2] = IntStream.range(0, double1.length).filter(i -> double1[i] > 0.0 && double2[i] > 0.0).count()
                / (double) double1.length;
        // Now we re-sort the counts and compute the top-N sims.  The default sort for counts is by highest-to-lowest score.
        Collections.sort(counts1);
        Collections.sort(counts2);
        int idx = 3;
        for (int nVal : this.nValueList) {
            if (counts1.size() >= nVal && counts2.size() >= nVal) {
                // Here both score lists have enough members to do a comparison.  Note that this also means both key
                // sets have the same number of members-- nVal.
                Set<String> keys1 = counts1.subList(0, nVal).stream().map(x -> x.getKey()).collect(Collectors.toSet());
                Set<String> keys2 = counts2.subList(0, nVal).stream().map(x -> x.getKey()).collect(Collectors.toSet());
                long common = keys1.stream().filter(x -> keys2.contains(x)).count();
                retVal[idx] = common / (double) nVal;
            }
            idx++;
        }
        return retVal;
    }


}
