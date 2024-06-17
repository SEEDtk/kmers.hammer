/**
 *
 */
package org.theseed.reports.eval;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.SortedMap;

import org.theseed.proteins.hammer.SummaryMap;

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
    protected final void summarizeFile(String name, SortedMap<String, SummaryMap> weightMaps) throws IOException {
        // Here we need to analyze the hits and produce an output line for each sample.  We process the
        // samples in order, and because the sample map is sorted, the samples will be presented in the same
        // order for each report.
        for (var sampEntry : weightMaps.entrySet()) {
            String sampleId = sampEntry.getKey();
            SummaryMap counters = sampEntry.getValue();
            SampleDescriptor desc = this.getSample(sampleId);
            if (desc == null)
                throw new IOException("Invalid sample ID " + sampleId + " found in report " + name + ".");
            // Get the expected best repgen.
            String expected = desc.getRepId();
            SummaryMap.Count counter0 = counters.findCounter(expected);
            double expectedCount = counter0.getCount();
            int expectedRoleCount = counter0.getRoleCount();
            // Count it as the best so far.
            String best = expected;
            double bestCount = expectedCount;
            // Accumulate bad hits here.
            double badCount = 0.0;
            int badRoleCount = 0;
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
                    int newBadRoleCount = counter.getRoleCount();
                    if (newBadRoleCount > badRoleCount)
                        badRoleCount = newBadRoleCount;
                }
            }
            this.processSampleData(name, desc, expectedCount, expectedRoleCount, best, bestCount, badCount, badRoleCount);
        }
    }

    /**
     * Process the summary analysis for a sample in a report.
     *
     * @param name				sample name
     * @param desc				sample descriptor
     * @param expectedCount		number of hits for expected repgen
     * @param expectedRoleCount number of roles for expected repgen
     * @param best				ID of best repgen
     * @param bestCount			number of hits for best repgen
     * @param badCount			number of hits for unexpected repgens
     * @param badRoleCount 		maximum role count for unexpected repgens
     */
    protected abstract void processSampleData(String name, SampleDescriptor desc, double expectedCount, int expectedRoleCount, String best,
            double bestCount, double badCount, int badRoleCount);

}
