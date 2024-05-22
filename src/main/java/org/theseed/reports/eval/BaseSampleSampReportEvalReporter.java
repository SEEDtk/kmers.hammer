/**
 *
 */
package org.theseed.reports.eval;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.SortedMap;

import org.theseed.stats.WeightMap;

/**
 * This report computes the top repgen hit results for each sample in a report, which gives us an indication of how each
 * sample performs in each report.  The subclass decides how to output the data.
 *
 * @author Bruce Parrello
 *
 */
public abstract class BaseSampleSampReportEvalReporter extends SummarySampReportEvalReporter {

    /**
     * Construct a sampling report for the specified controlling command processor and the
     * specified report writer.
     *
     * @param processor		controlling command processor
     * @param writer		output writer for report
     */
    public BaseSampleSampReportEvalReporter(IParms processor, PrintWriter writer) {
        super(processor, writer);
    }

    @Override
    protected final void summarizeFile(String name, SortedMap<String, WeightMap> weightMaps) throws IOException {
        // Here we need to analyze the hits and produce an output line for each sample.  We process the
        // samples in order, and because the sample map is sorted, the samples will be presented in the same
        // order for each report.
        for (var sampEntry : weightMaps.entrySet()) {
            String sampleId = sampEntry.getKey();
            WeightMap counters = sampEntry.getValue();
            SampleDescriptor desc = this.getSample(sampleId);
            if (desc == null)
                throw new IOException("Invalid sample ID " + sampleId + " found in report " + name + ".");
            // Get the expected best repgen.
            String expected = desc.getRepId();
            double expectedCount = counters.getCount(expected);
            // Count it as the best so far.
            String best = expected;
            double bestCount = expectedCount;
            // Accumulate bad hits here.
            double badCount = 0.0;
            // Loop through the counts.
            for (var counter : counters.counts()) {
                // Get the repgen ID and count.
                String other = counter.getKey();
                double count = counter.getCount();
                if (! other.equals(expected)) {
                    // Here we have bad hits.
                    if (count > bestCount) {
                        best = other;
                        bestCount = count;
                    }
                    badCount += count;
                }
            }
            this.processSampleData(name, desc, expectedCount, best, bestCount, badCount);
        }
    }

    /**
     * Process the summary analysis for a sample in a report.
     *
     * @param name				sample name
     * @param desc				sample descriptor
     * @param expectedCount		number of hits for expected repgen
     * @param best				ID of best repgen
     * @param bestCount			number of hits for best repgen
     * @param badCount			number of hits for unexpected repgens
     */
    protected abstract void processSampleData(String name, SampleDescriptor desc, double expectedCount, String best,
            double bestCount, double badCount);

}
