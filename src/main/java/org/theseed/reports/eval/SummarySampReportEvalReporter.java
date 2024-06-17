/**
 *
 */
package org.theseed.reports.eval;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.HashMap;

import org.theseed.proteins.hammer.SummaryMap;

/**
 * This reporter tabulates the good and bad hits for each sample and provides hooks for the subclass to output
 * the resulting analysis.
 *
 * @author Bruce Parrello
 *
 */
public abstract class SummarySampReportEvalReporter extends SampReportEvalReporter {

    // FIELDS
    /** name of the current report */
    private String reportName;
    /** map of samples to weight maps for the current report */
    private SortedMap<String, SummaryMap> countMap;
    /** map of repgen IDs to names */
    private Map<String, String> repMap;

    /**
     * Construct a quality report.
     *
     * @param processor		controlling command processor
     * @param writer		output print writer
     */
    public SummarySampReportEvalReporter(IParms processor, PrintWriter writer) {
        super(processor, writer);
        // Create the maps.
        this.countMap = new TreeMap<String, SummaryMap>();
        this.repMap = new HashMap<String, String>();
        // Denote no report is active.
        this.reportName = null;
    }

    @Override
    public final void openFile(String name) {
        this.reportName = name;
        // Erase all the current counts.
        this.countMap.clear();
        // Perform any required file initialization.
        this.initFile(name);
    }

    /**
     * Initialize for a new report file.
     *
     * @param name		new file name
     */
    protected abstract void initFile(String name);

    @Override
    public final void recordHits(SampleDescriptor desc, String repId, String repName, double count, int roleCount) {
        String sampleId = desc.getSampleId();
        // Get the weight map for this sample and update it.
        SummaryMap counters = this.countMap.computeIfAbsent(sampleId, x -> new SummaryMap());
        counters.count(repId, count, roleCount);
        // Insure we have the name saved for the repgen ID.
        this.repMap.computeIfAbsent(repId, x -> repName);
    }

    @Override
    public final void closeFile() throws IOException {
        this.summarizeFile(this.reportName, this.countMap);
    }

    /**
     * Summarize and output the results for the current report.
     *
     * @param name			current report
     * @param countMap2		map of sample names to repgen weight totals
     *
     * @throws IOException
     */
    protected abstract void summarizeFile(String name, SortedMap<String, SummaryMap> countMap2) throws IOException;

    /**
     * @return the name of a repgen genome
     *
     * @param repId		ID of the repgen genome
     */
    public String getRepName(String repId) {
        return this.repMap.getOrDefault(repId, "<unknown genome>");
    }

}
