/**
 *
 */
package org.theseed.reports.eval;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.SortedMap;

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
public class QualitySampReportEvalReporter extends SummarySampReportEvalReporter {

    /**
     * Construct a quality report.
     *
     * @param processor		controlling command processor
     * @param writer		output print writer
     */
    public QualitySampReportEvalReporter(IParms processor, PrintWriter writer) {
        super(processor, writer);
        // Write the header line.
        writer.println("report_name\tsamples\tgood_samples\tbad_samples\tgood_hits\tgood_hit%\tbad_hits\tbad_hit%\tmean_bad%");
    }

    @Override
    protected void initFile(String name) {
    }

    @Override
    protected void summarizeFile(String reportName, SortedMap<String, WeightMap> countMap) throws IOException {
        // Total up the data for this file.
        double totalGood = 0.0;
        double totalBad = 0.0;
        int goodSamples = 0;
        int badSamples = 0;
        double totalPct = 0.0;
        // We loop through the samples.  For each sample, we process its weight map and determine the
        // highest-weight repgen.  If it's the expected one for the sample, the sample is good.
        for (var sampEntry : countMap.entrySet()) {
            String sampleId = sampEntry.getKey();
            WeightMap counters = sampEntry.getValue();
            // Get the expected repgen.
            SampleDescriptor desc = this.getSample(sampleId);
            if (desc == null)
                throw new IOException("Invalid sample ID " + sampleId + " found in file " + reportName + ".");
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
            if (bad > 0.0)
                totalPct += (bad * 100.0) / (good + bad);
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
        double meanBadPct = 0.0;
        if (totSamples > 0)
            meanBadPct = totalPct / totSamples;
        writer.println(reportName + "\t" + totSamples + "\t" + goodSamples + "\t" + badSamples
                + "\t" + totalGood + "\t" + pctGood + "\t" + totalBad + "\t" + pctBad + "\t" + meanBadPct);
    }

    @Override
    public void closeReport() {
    }

}
