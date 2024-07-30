/**
 *
 */
package org.theseed.reports.eval;

import java.io.PrintWriter;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.commons.lang3.StringUtils;
import org.theseed.proteins.hammer.SummaryMap.Count;

/**
 * This report lists each repgen hit by a sample in a report, along with the distance to the repgen (if any),
 * and the hit counts and percentage.
 *
 * @author Bruce Parrello
 *
 */
public class DetailSampReportEvalReporter extends DetailBaseSampReportEvalReporter {

    public DetailSampReportEvalReporter(IParms processor, PrintWriter writer) {
        super(processor, writer);
    }

    @Override
    protected void writeHeader(List<String> roles) {
        writer.println("file_name\tsample_id\tsample_genome\trepgen_id\trepgen_name\thits\tpct_hits\trole_count\tdistance\t"
                + StringUtils.join(roles, "\t"));
    }

    @Override
    protected void outputDetail(String reportName, String sampleId, String sampleGenomeId, double count,
            double pctCount, Count repHit) {
        // Get the representative genome data.
        String repId = repHit.getKey();
        int roleCount = repHit.getNumRoles();
        String repName = this.getRepName(repId);
        // Compute the distance.
        double distance = this.getDistance(sampleGenomeId, repId);
        // Get the role counts.
        List<String> roleList = this.getRoles();
        String roleCountLine = roleList.stream().map(x -> Double.toString(repHit.getRoleCount(x))).collect(Collectors.joining("\t"));
        this.println(reportName + "\t" + sampleId + "\t" + sampleGenomeId + "\t" + repId + "\t"
                + repName + "\t" + count + "\t" + pctCount + "\t" + roleCount + "\t" + distance
                + "\t" + roleCountLine);
    }

    @Override
    public void closeReport() {
    }

}
