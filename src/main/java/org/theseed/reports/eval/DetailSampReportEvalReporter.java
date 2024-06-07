/**
 *
 */
package org.theseed.reports.eval;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

import org.theseed.stats.WeightMap;

/**
 * This report lists each repgen hit by a sample in a report, along with the distance to the repgen (if any),
 * and the hit counts and percentage.
 *
 * @author Bruce Parrello
 *
 */
public class DetailSampReportEvalReporter extends SampReportEvalReporter {

    // FIELDS
    /** name of the current file */
    private String fileName;
    /** map of sample IDs to sample descriptors for current file */
    private Map<String, SampleDescriptor> sampleMap;
    /** map of repgen IDs to names for the current file */
    private Map<String, String> repNameMap;
    /** number of hits for each repgen ID in the current file, keyed by sample ID */
    private Map<String, WeightMap> repHitMap;

    /**
     * Construct a detail sample evaluation report.
     *
     * @param processor		controlling command processor
     * @param writer		output print writer
     */
    public DetailSampReportEvalReporter(IParms processor, PrintWriter writer) {
        super(processor, writer);
        this.sampleMap = new HashMap<String, SampleDescriptor>();
        this.repHitMap = new TreeMap<String, WeightMap>();
        this.repNameMap = new HashMap<String, String>();
        // Write out the report header.
        writer.println("file_name\tsample_id\tsample_genome\trepgen_id\trepgen_name\thits\tpct_hits\tdistance");
    }

    @Override
    public void openFile(String name) {
        // Initialize all the data structures for the file.
        this.fileName = name;
        this.sampleMap.clear();
        this.repHitMap.clear();
        this.repNameMap.clear();
    }

    @Override
    public void recordHits(SampleDescriptor desc, String repId, String repName, double count) throws IOException {
        // Store the basic data.
        this.repNameMap.put(repId, repName);
        final String sampleId = desc.getSampleId();
        this.sampleMap.put(sampleId, desc);
        // Record the hits.
        WeightMap repHits = this.repHitMap.computeIfAbsent(sampleId, x -> new WeightMap());
        repHits.count(repId, count);
    }

    @Override
    public void closeFile() throws IOException {
        // Here we have processed a whole file.  Now we need to output the data.  We loop in
        // two levels, starting with the sample IDs.
        for (var sampleHitEntry : this.repHitMap.entrySet()) {
            // Get the sample data.
            String sampleId = sampleHitEntry.getKey();
            SampleDescriptor sample = this.sampleMap.get(sampleId);
            String sampleGenomeId = sample.getGenomeId();
            // Now get the map of hits per repgen for this sample.
            WeightMap repHits = sampleHitEntry.getValue();
            // Compute the total hits for this sample.
            double total = repHits.sum();
            // Loop through the hits in sorted order.
            for (var repHit : repHits.sortedCounts()) {
                String repId = repHit.getKey();
                double count = repHit.getCount();
                String repName = this.repNameMap.get(repId);
                // Compute the hit percentage.
                double pctCount = 0.0;
                if (count > 0.0)
                    pctCount = count * 100.0 / total;
                // Compute the distance.
                double distance = this.getDistance(sampleGenomeId, repId);
                // Write the report line.
                this.println(this.fileName + "\t" + sampleId + "\t" + sampleGenomeId + "\t" + repId + "\t"
                        + repName + "\t" + count + "\t" + pctCount + "\t" + distance);
            }
        }
    }

    @Override
    public void closeReport() {
    }

}
