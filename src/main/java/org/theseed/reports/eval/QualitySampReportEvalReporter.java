/**
 *
 */
package org.theseed.reports.eval;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.Map;
import java.util.HashMap;

import org.theseed.stats.WeightMap;

/**
 * This reporter tabulates the good and bad hits for each sample and outputs a comparison of the various bin-reports.
 * For each bin-report, it outputs the number of good hits, the number of bad hits, the percentages for each, the number
 * of correct samples, and the number of incorrect samples.  In this case, a sample is correct if the good hits outnumber
 * the bad ones.
 *
 * @author Bruce Parrello
 *
 */
public class QualitySampReportEvalReporter extends SampReportEvalReporter {

    // FIELDS
    /** name of the current report */
    private String reportName;
    /** map of samples to weight maps for the current report */
    private Map<String, WeightMap> countMap;

    /**
     * Construct a quality report.
     *
     * @param processor		controlling command processor
     * @param writer		output print writer
     */
    public QualitySampReportEvalReporter(IParms processor, PrintWriter writer) {
        super(processor, writer);
        // Create the count map.
        this.countMap = new HashMap<String, WeightMap>();
        // Denote no report is active.
        this.reportName = null;
        // Write the header line.
        writer.println("report_name\tsamples\tgood_samples\tbad_samples\tgood_hits\tgood_hit%\tbad_hits\tbad_hit%");
    }

    @Override
    public void openFile(String name) {
        this.reportName = name;
        // Erase all the current counts.
        this.countMap.clear();
    }

    @Override
    public void recordHits(SampleDescriptor desc, String repId, String repName, double count) {
        String sampleId = desc.getSampleId();
        // Get the weight map for this sample.
        WeightMap counters = this.countMap.computeIfAbsent(sampleId, x -> new WeightMap());
        counters.count(repId, count);
    }

    @Override
    public void closeFile() throws IOException {
        // Total up the data for this file.
        double totalGood = 0.0;
        double totalBad = 0.0;
        int goodSamples = 0;
        int badSamples = 0;
        // We loop through the samples.  For each sample, we process its weight map and determine the
        // highest-weight repgen.  If it's the expected one for the sample, the sample is good.
        for (var sampEntry : this.countMap.entrySet()) {
            String sampleId = sampEntry.getKey();
            WeightMap counters = sampEntry.getValue();
            // Get the expected repgen.
            SampleDescriptor desc = this.getSample(sampleId);
            if (desc == null)
                throw new IOException("Invalid sample ID " + sampleId + " found in file " + this.reportName + ".");
            String expected = desc.getRepId();
            // Get the number of hits on the expected repgen and set up to count.
            double good = counters.getCount(expected);
            double bad = 0.0;
            boolean goodSample = true;
            // Now loop through the other counts.
            for (var counter : counters.counts()) {
                String repId = counter.getKey();
                if (! repId.equals(expected)) {
                    double badCount = counter.getCount();
                    // If this incorrect result is not less than the expected result, this is a bad sample.
                    if (badCount >= good)
                        goodSample = false;
                    bad += badCount;
                }
            }
            // Accumulate this sample's results.
            if (goodSample)
                goodSamples++;
            else
                badSamples++;
            totalGood += good;
            totalBad += bad;
        }
        // Compute the percents and totals.
        int totSamples = goodSamples + badSamples;
        double totHits = totalGood + totalBad;
        double pctGood = 0.0;
        double pctBad = 0.0;
        if (totHits > 0.0) {
            pctGood = totalGood * 100.0 / totHits;
            pctBad = totalBad * 100.0 / totHits;
        }
        writer.println(this.reportName + "\t" + totSamples + "\t" + goodSamples + "\t" + badSamples
                + "\t" + totalGood + "\t" + pctGood + "\t" + totalBad + "\t" + pctBad);
    }

    @Override
    public void closeReport() {
    }


}
