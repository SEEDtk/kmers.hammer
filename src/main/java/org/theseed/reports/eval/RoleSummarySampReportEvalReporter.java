/**
 *
 */
package org.theseed.reports.eval;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedMap;
import java.util.stream.Collectors;

import org.apache.commons.lang3.StringUtils;
import org.theseed.proteins.hammer.SummaryMap;

/**
 * This report is designed to determine the performance of different roles in the hammer analysis of a sample.  For each
 * sample in a report file, it will output the good and bad hits by role.
 *
 * @author Bruce Parrello
 *
 */
public class RoleSummarySampReportEvalReporter extends SummarySampReportEvalReporter {

    // FIELDS
    /** list of roles of interest */
    private List<String> roleList;

    /**
     * Construct a role-based summary report.
     *
     * @param processor		controlling command processor
     * @param writer		output print writer
     */
    public RoleSummarySampReportEvalReporter(IParms processor, PrintWriter writer) {
        super(processor, writer);
        // Get the role set.
        this.roleList = new ArrayList<String>(processor.getRoleSet());
        // Write the header.
        writer.println("file_name\tsample_id\ttype\tcount\t" + StringUtils.join(this.roleList, "\t"));
    }

    @Override
    protected void initFile(String name) {
    }

    @Override
    protected void summarizeFile(String name, SortedMap<String, SummaryMap> countMap) throws IOException {
        // We keep some totals in here.
        SummaryMap totalCounts = new SummaryMap();
        for (var sampEntry : countMap.entrySet()) {
            String sampleId = sampEntry.getKey();
            SummaryMap counters = sampEntry.getValue();
            // Get the expected repgen.
            SampleDescriptor desc = this.getSample(sampleId);
            if (desc == null)
                throw new IOException("Invalid sample ID " + sampleId + " found in file " + name + ".");
            String expected = desc.getRepId();
            // Now we loop through the repgens for this sample, accumulating the good and bad data.
            SummaryMap summaryCounts = new SummaryMap();
            for (var repCount : counters.counts()) {
                String repGenId = repCount.getKey();
                String type = (repGenId.equals(expected) ? "good" : "bad");
                summaryCounts.merge(type, repCount);
                totalCounts.merge(type, repCount);
            }
            // Finally, we write out the good and bad data.
            this.writeSummaryCounts(name, sampleId, summaryCounts);
        }
        // Now write out the totals for the file.
        this.writeSummaryCounts(name, "TOTALS", totalCounts);

    }

    /**
     * Write the good and bad summary counts from a summary map.
     *
     * @param fileName		name of the current report file
     * @param rowName		row label (usually the sample ID)
     * @param summaryMap	summary map to write
     */
    private void writeSummaryCounts(String fileName, String rowName, SummaryMap summaryMap) {
        var goodCount = summaryMap.findCounter("good");
        if (goodCount != null)
            this.writeLine(fileName, rowName, goodCount);
        var badCount = summaryMap.findCounter("bad");
        if (badCount != null)
            this.writeLine(fileName, rowName, badCount);
    }

    /**
     * Write a line of output for this report.
     *
     * @param name			current file name
     * @param sampleId		relevant sample ID
     * @param counter		counter containing the data
     */
    private void writeLine(String name, String sampleId, SummaryMap.Count counter) {
        // Construct the role counts.
        String roleCountLine = this.roleList.stream().map(x -> Double.toString(counter.getRoleCount(x))).collect(Collectors.joining("\t"));
        writer.println(name + "\t" + sampleId + "\t" + counter.getKey() + "\t" + counter.getCount() + "\t" + roleCountLine);
    }

    @Override
    public void closeReport() {
    }

}
