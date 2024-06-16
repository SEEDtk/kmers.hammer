/**
 *
 */
package org.theseed.reports.eval;

import java.io.PrintWriter;

/**
 * This report outputs the top repgen hit results for each sample in a report, which gives us an indication of how each
 * sample performs in each report.
 *
 * @author Bruce Parrello
 *
 */
public class SamplingSampReportEvalReporter extends BaseSampleSampReportEvalReporter {

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
        writer.println("report_name\tsample_id\tsample_genome\texpected_rep\texpected_name\texpected_hits\texpected_hit%\tbest_rep\tbest_name\tbest_hits\tbest_hit%\tbad_hits\tbad_hits%\tdistance");
    }

    @Override
    protected void initFile(String name) {
    }

    @Override
    protected void processSampleData(String name, SampleDescriptor desc, double expectedCount, String best,
            double bestCount, double badCount) {
        // Get the data for the sample.
        String sampleId = desc.getSampleId();
        String sampleGenome = desc.getGenomeId();
        // We need to compute the percents for this sample.
        double totalCount = expectedCount + badCount;
        double expectedPct = 0.0;
        double badPct = 0.0;
        double bestPct = 0.0;
        if (totalCount > 0.0) {
            expectedPct = expectedCount * 100.0 / totalCount;
            badPct = badCount * 100.0 / totalCount;
            bestPct = bestCount * 100.0 / totalCount;
        }
        // Get the data for the expected repgen.
        String expected = desc.getRepId();
        String expectedName = this.getRepName(expected);
        // Get the distance to the expected repgen.
        double distance = this.getDistance(sampleGenome, expected);
        // Finally, get the name of the best genome.
        String bestName = this.getRepName(best);
        // Write the sample data.
        writer.println(name + "\t" + sampleId + "\t" + sampleGenome + "\t" + expected + "\t" + expectedName + "\t" + expectedCount + "\t" + expectedPct
                + "\t" + best + "\t" + bestName + "\t" + bestCount + "\t" + bestPct + "\t" + badCount + "\t" + badPct + "\t" + distance);
    }

    @Override
    public void closeReport() {
    }

}
