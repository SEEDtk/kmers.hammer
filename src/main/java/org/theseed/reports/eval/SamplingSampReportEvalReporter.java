/**
 *
 */
package org.theseed.reports.eval;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.SortedMap;

import org.theseed.stats.WeightMap;

/**
 * This report outputs the top repgen hit results for each sample in a report, which gives us an indication of how each
 * sample performs in each report.
 *
 * @author Bruce Parrello
 *
 */
public class SamplingSampReportEvalReporter extends SummarySampReportEvalReporter {

    /**
     * Construct a sampling report for the specified controlling command processor and the
     * specified report writer.
     *
     * @param processor		controlling command processor
     * @param writer		output writer for report
     */
    public SamplingSampReportEvalReporter(IParms processor, PrintWriter writer) {
        super(processor, writer);
        // Write the report header.
        writer.println("report_name\tsample_id\tsample_genome\texpected_rep\texpected_name\texpected_hits\texpected_hit%\tbest_rep\tbest_name\tbest_hits\tbest_hit%\tbad_hits\tbad_hits%");
    }

    @Override
    protected void initFile(String name) {
    }

    @Override
    protected void summarizeFile(String name, SortedMap<String, WeightMap> weightMaps) throws IOException {
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
            String expectedName = this.getRepName(expected);
            double expectedCount = counters.getCount(expected);
            // Get the sample's genome ID.
            String sampleGenome = desc.getGenomeId();
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
            // Now we need to compute the percents for this sample.
            double totalCount = expectedCount + badCount;
            double expectedPct = 0.0;
            double badPct = 0.0;
            double bestPct = 0.0;
            if (totalCount > 0.0) {
                expectedPct = expectedCount * 100.0 / totalCount;
                badPct = badCount * 100.0 / totalCount;
                bestPct = bestCount * 100.0 / totalCount;
            }
            // Finally, get the name of the best genome.
            String bestName = this.getRepName(best);
            // Write the sample data.
            writer.println(name + "\t" + sampleId + "\t" + sampleGenome + "\t" + expected + "\t" + expectedName + "\t" + expectedCount + "\t" + expectedPct
                    + "\t" + best + "\t" + bestName + "\t" + bestCount + "\t" + bestPct + "\t" + badCount + "\t" + badPct);
        }
    }

    @Override
    public void closeReport() {
    }

}
