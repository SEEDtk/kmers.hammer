/**
 *
 */
package org.theseed.reports.eval;

import java.io.PrintWriter;
import java.util.SortedMap;
import java.util.TreeMap;

/**
 * This report lists the best bin report for each sample.
 *
 * @author Bruce Parrello
 *
 */
public class RatingSampReportEvalReporter extends BaseSampleSampReportEvalReporter {

    // FIELDS
    /** map of sample IDs to best report for each sample */
    private SortedMap<String, BestReport> bestMap;

    /**
     * This class tracks the best report for each sample.  It is sorted by bad-hit-percent (lower sorts first),
     * then correct repgen (correct sorts before incorrect), then the total number of bad hits (smaller is better).
     *
     * @author Bruce Parrello
     *
     */
    protected static class BestReport implements Comparable<BestReport> {

        /** name of the best report */
        private String reportName;
        /** percent of bad hits */
        private double badPercent;
        /** TRUE if the best repgen is the expected repgen */
        private boolean correctRepGen;
        /** number of bad hits */
        private double badCount;
        /** total number of hits */
        private double hitCount;

        /**
         * Construct a new best-report object from the report data for a sample.
         *
         * @param name		name of report
         * @param badPct	percent of hits that were bad
         * @param badCount	number of hits that were bad
         * @param hitCount	total number of hits
         * @param expected	expected repgen ID
         * @param best		actual repgen ID
         */
        protected BestReport(String name, double badPct, double badCount, double hitCount, String expected, String best) {
            this.reportName = name;
            this.badPercent = badPct;
            this.badCount = badCount;
            this.correctRepGen = expected.equals(best);
            this.hitCount = hitCount;
        }

        @Override
        public int compareTo(BestReport o) {
            int retVal = Double.compare(this.badPercent, o.badPercent);
            if (retVal == 0) {
                retVal = (this.correctRepGen ? (o.correctRepGen ? 0 : -1) : (o.correctRepGen ? 1 : 0));
                if (retVal == 0) {
                    retVal = Double.compare(this.badCount, o.badCount);
                    if (retVal == 0)
                        retVal = this.reportName.compareTo(o.reportName);
                }
            }
            return retVal;
        }

        /**
         * @return the report name
         */
        public String getReportName() {
            return this.reportName;
        }

        /**
         * @return the percent of hits that were bad
         */
        public double getBadPercent() {
            return this.badPercent;
        }

        /**
         * @return TRUE if the best repgen was the expected one, else FALSE
         */
        public boolean isCorrectRepGen() {
            return this.correctRepGen;
        }

        /**
         * @return the number of hits that were bad
         */
        public double getBadCount() {
            return this.badCount;
        }

        /**
         * @return the total number of hits
         */
        public double getHitCount() {
            return this.hitCount;
        }
    }

    /**
     * Construct a new rating sample reporter.
     *
     * @param processor		controlling command processor
     * @param writer		output report writer
     */
    public RatingSampReportEvalReporter(IParms processor, PrintWriter writer) {
        super(processor, writer);
        // Create the tracking map.
        this.bestMap = new TreeMap<String, BestReport>();
    }

    @Override
    protected void processSampleData(String name, SampleDescriptor desc, double expectedCount, int expectedRoleCount, String best,
            double bestCount, int bestRoleCount, double badCount, int badRoleCount) {
        // Compute the bad-hit percentage.
        double badPercent = 0.0;
        if (badCount > 0)
            badPercent = badCount * 100.0 / (badCount + expectedCount);
        // Compute the best-report tracker.
        BestReport reportData = new BestReport(name, badPercent, badCount, badCount + expectedCount, desc.getRepId(), best);
        // Compare it to the existing value.  Store it if it is better than the old value.
        String sampleId = desc.getSampleId();
        BestReport oldData = this.bestMap.get(sampleId);
        if (oldData == null || reportData.compareTo(oldData) < 0)
            this.bestMap.put(sampleId, reportData);
    }

    @Override
    protected void initFile(String name) {
    }

    @Override
    public void closeReport() {
        // All done, so write the results.  We start with the header.
        writer.println("sample_id\tbest_report\thit_count\tbad_hits\tbad_percent\tcorrect_repgen");
        // Now, loop through the best-report map.
        for (var sampEntry : this.bestMap.entrySet()) {
            String sampleId = sampEntry.getKey();
            BestReport reportData = sampEntry.getValue();
            writer.println(sampleId + "\t" + reportData.getReportName() + "\t" + reportData.getHitCount() + "\t" +  reportData.getBadCount()
                    + "\t" + reportData.getBadPercent() + "\t" + (reportData.isCorrectRepGen() ? "Y" : ""));
        }
    }
}
