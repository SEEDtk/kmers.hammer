/**
 *
 */
package org.theseed.reports.eval;

import java.io.PrintWriter;

/**
 * This report lists the best bin report for each sample.
 *
 * @author Bruce Parrello
 *
 */
public class RatingSampReportEvalReporter extends BaseSampleSampReportEvalReporter {

    /**
     * @param processor
     * @param writer
     */
    public RatingSampReportEvalReporter(IParms processor, PrintWriter writer) {
        super(processor, writer);
        // TODO Auto-generated constructor stub
    }

    @Override
    protected void processSampleData(String name, SampleDescriptor desc, double expectedCount, String best,
            double bestCount, double badCount) {
        // TODO code for processSampleData

    }

    @Override
    protected void initFile(String name) {
        // TODO code for initFile

    }

    @Override
    public void closeReport() {
        // TODO code for closeReport

    }
    // FIELDS
    // TODO data members for RatingSampReportEvalReporter

    // TODO constructors and methods for RatingSampReportEvalReporter
}
