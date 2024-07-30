/**
 *
 */
package org.theseed.reports.eval;

import java.io.PrintWriter;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.proteins.hammer.SummaryMap.Count;

/**
 * This report organizes the detail data by sample/repgen pair and outputs the results showing the
 * differences between the hit percentages for each report.  The output is sorted by distance, which
 * allows for loading the result into a spreadsheet for graphing.
 *
 * @author Bruce Parrello
 *
 */
public class CompareSampReportEvalReporter extends DetailBaseSampReportEvalReporter {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(CompareSampReportEvalReporter.class);
    /** hit summary map; sample -> genome -> summary object */
    private Map<String, Map<String, HitSummary>> summaryMap;
    /** list of report names */
    private SortedSet<String> reportNames;

    /**
     * This class contains the output data for each repgen/sample pairing.
     */
    protected class HitSummary implements Comparable<HitSummary> {

        /** sample ID */
        private String sampleId;
        /** sample genome ID */
        private String sampleGenomeId;
        /** repgen ID */
        private String repId;
        /** distance */
        private double distance;
        /** map of report names to hit percents */
        private TreeMap<String, Double> pctMap;

        /**
         * Create a new hit-summary object.
         *
         * @param sample	ID of the relevant sample
         * @param genome	Id of the sample's genome
         * @param rep		ID of the relevant representative genome
         * @param dist		distance from the sample genome to the representative genome
         */
        protected HitSummary(String sample, String genome, String rep, double dist) {
            this.sampleId = sample;
            this.sampleGenomeId = genome;
            this.repId = rep;
            this.distance = dist;
            this.pctMap = new TreeMap<String, Double>();
        }

        /**
         * Stored a hit percentage in this hit summary.
         *
         * @param reportName	name of relevant report
         * @param percent		percentage to store
         */
        protected void update(String reportName, double percent) {
            this.pctMap.put(reportName, percent);
        }

        @Override
        public int compareTo(HitSummary o) {
            // Sort ascending by distance.
            int retVal = Double.compare(this.distance, o.distance);
            if (retVal == 0) {
                // Break ties on sample ID, then repgen ID.
                retVal = this.sampleId.compareTo(o.sampleId);
                if (retVal == 0)
                    retVal = this.repId.compareTo(o.repId);
            }
            return retVal;
        }

        /**
         * @return the sample ID
         */
        public String getSampleId() {
            return this.sampleId;
        }

        /**
         * @return the representative genome ID
         */
        public String getRepId() {
            return this.repId;
        }

        /**
         * @return the distance between the sample genome and the representative genome
         */
        public double getDistance() {
            return this.distance;
        }

        /**
         * @return the hit percentage for the specified report
         *
         * @param reportName	name of the report of interest
         */
        public double getHitPercent(String reportName) {
            Double retVal = this.pctMap.get(reportName);
            // No hits means zero percent.
            if (retVal == null)
                retVal = 0.0;
            return retVal;
        }

        /**
         * @return the ID of the genome represented by the sample
         */
        public String getSampleGenomeId() {
            return this.sampleGenomeId;
        }

    }

    /**
     * Construct a comparison report.
     *
     * @param processor		controlling command processor
     * @param writer		output print writer
     */
    public CompareSampReportEvalReporter(IParms processor, PrintWriter writer) {
        super(processor, writer);
        this.reportNames = new TreeSet<String>();
        this.summaryMap = new TreeMap<String, Map<String, HitSummary>>();
    }

    @Override
    protected void writeHeader(List<String> roles) {
        // We can't write the header until the end.
    }

    @Override
    protected void outputDetail(String reportName, String sampleId, String sampleGenomeId, double count,
            double pctCount, Count repHit) {
        // First, insure we remember this report name.
        this.reportNames.add(reportName);
        // Get the map for this sample ID.
        Map<String, HitSummary> sampleMap = this.summaryMap.computeIfAbsent(sampleId, x -> new HashMap<String, HitSummary>());
        // Get the summary for this representative genome.
        String repId = repHit.getKey();
        double distance = this.getDistance(sampleGenomeId, repId);
        HitSummary hitSummary = sampleMap.computeIfAbsent(repId, x -> new HitSummary(sampleId, sampleGenomeId, repId, distance));
        // Add this report's data to the hit summary.
        hitSummary.update(reportName, pctCount);
    }

    @Override
    public void closeReport() {
        // Form a header line from the saved report names.
        String reportHeader = StringUtils.join(this.reportNames, '\t');
        writer.println("sample_id\tsample_genome\trep_id\trep_name\tdistance\t" + reportHeader);
        // Loop through the hit summaries, writing the report.  Here we collect them together and do the sort.
        List<HitSummary> summaries = this.summaryMap.values().stream().flatMap(x -> x.values().stream())
                .sorted().collect(Collectors.toList());
        // We will build the output lines in here.
        StringBuilder buffer = new StringBuilder(100);
        for (HitSummary summary : summaries) {
            buffer.setLength(0);
            buffer.append(summary.getSampleId());
            buffer.append('\t').append(summary.getSampleGenomeId());
            String repId = summary.getRepId();
            String repName = this.getRepName(repId);
            buffer.append('\t').append(repId).append('\t').append(repName);
            buffer.append('\t').append(summary.getDistance());
            for (String reportName : this.reportNames) {
                double hitPct = summary.getHitPercent(reportName);
                buffer.append('\t').append(hitPct);
            }
            writer.println(buffer.toString());
        }
    }

}
