/**
 *
 */
package org.theseed.reports.eval;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import org.theseed.counters.CountMap;
import org.theseed.proteins.hammer.SummaryMap;
import org.theseed.proteins.hammer.SummaryMap.Count;

/**
 * This base class is for reports that want to focus on detailed sample/repgen hit relationships.  It provides the hit percentage
 * and the raw counts.
 *
 * @author Bruce Parrello
 *
 */
public abstract class DetailBaseSampReportEvalReporter extends SampReportEvalReporter {

    // FIELDS
    /** name of the current file */
    private String fileName;
    /** map of sample IDs to sample descriptors for current file */
    private Map<String, SampleDescriptor> sampleMap;
    /** map of repgen IDs to names for the current file */
    private Map<String, String> repNameMap;
    /** number of hits for each repgen ID in the current file, keyed by sample ID */
    private Map<String, SummaryMap> repHitMap;
    /** list of roles of interest */
    private List<String> roleList;

    /**
     * Construct a detail sample evaluation report.
     *
     * @param processor		controlling command processor
     * @param writer		output print writer
     */
    public DetailBaseSampReportEvalReporter(IParms processor, PrintWriter writer) {
        super(processor, writer);
        this.sampleMap = new HashMap<String, SampleDescriptor>();
        this.repHitMap = new TreeMap<String, SummaryMap>();
        this.repNameMap = new HashMap<String, String>();
        this.roleList = new ArrayList<String>(processor.getRoleSet());
        // Write out the report header.
        this.writeHeader(this.roleList);
    }

    /**
     * Write out the header line for this report.
     *
     * @param roles		list of roles used
     */
    protected abstract void writeHeader(List<String> roles);

    @Override
    public void openFile(String name) {
        // Initialize all the data structures for the file.
        this.fileName = name;
        this.sampleMap.clear();
        this.repHitMap.clear();
    }

    @Override
    public void recordHits(SampleDescriptor desc, String repId, String repName, double count, int roleCount, CountMap<String> roleCounts) throws IOException {
        // Store the basic data.
        this.repNameMap.put(repId, repName);
        final String sampleId = desc.getSampleId();
        this.sampleMap.put(sampleId, desc);
        // Record the hits.
        SummaryMap repHits = this.repHitMap.computeIfAbsent(sampleId, x -> new SummaryMap());
        repHits.count(repId, count, roleCount, roleCounts);
    }

    @Override
    public void closeFile() throws IOException {
        // Here we have processed a whole file.  Now we need to output the data.
        // We loop in two levels, starting with the sample IDs.
        for (var sampleHitEntry : this.repHitMap.entrySet()) {
            // Get the sample data.
            String sampleId = sampleHitEntry.getKey();
            SampleDescriptor sample = this.sampleMap.get(sampleId);
            String sampleGenomeId = sample.getGenomeId();
            // Now get the map of hits per repgen for this sample.
            SummaryMap repHits = sampleHitEntry.getValue();
            // Compute the total hits for this sample.
            double total = repHits.sum();
            // Loop through the hits in sorted order.
            for (var repHit : repHits.sortedCounts()) {
                double count = repHit.getCount();
                // Compute the hit percentage.
                double pctCount = 0.0;
                if (count > 0.0)
                    pctCount = count * 100.0 / total;
                // Write the report line.
                this.outputDetail(this.fileName, sampleId, sampleGenomeId, count, pctCount, repHit);
            }
        }
    }

    /**
     * Output the details for a specific hit relationship between a sample and a representative genome.
     *
     * @param reportName		name of the relevant report
     * @param sampleId			ID of the sample
     * @param sampleGenomeId	ID of the sample's genome
     * @param count				total strength of the hits
     * @param pctCount			percent of all the hits for the sample belonging to this repgen
     * @param repHit			count object for the hits
     */
    protected abstract void outputDetail(String reportName, String sampleId, String sampleGenomeId, double count,
            double pctCount, Count repHit);

    /**
     * @return the name of a representative genome
     *
     * @param repId		ID of the genome of interest
     */
    public String getRepName(String repId) {
        String retVal = this.repNameMap.get(repId);
        if (retVal == null)
            retVal = "<unknown repgen " + repId + ">";
        return retVal;
    }

    /**
     * @return the list of roles
     */
    public List<String> getRoles() {
        return this.roleList;
    }

}
