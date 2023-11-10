/**
 *
 */
package org.theseed.reports;

import java.io.IOException;
import java.util.Map;
import java.util.stream.Collectors;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.proteins.hammer.HammerDb.Source;

/**
 * This report produces a simple listing of the hammers along with the role associated
 * with each one.
 *
 * @author Bruce Parrello
 *
 */
public class GenomeHammerReport extends HammerReport {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(GenomeHammerReport.class);


    public GenomeHammerReport(IParms processor) throws ParseFailureException, IOException {
        super(processor);
    }

    @Override
    protected void initReport() {
        // Write the header line.
        this.printFields("hammer", "fid", "rep_id", "rep_name", "role_id", "role");
    }

    @Override
    public void processGenome(String repId, String repName, Map<String, Source> hammerMap) {
        // Get the source representative genome.
        Genome repGenome = this.getRepGenome(repId);
        // Loop through the hammers, producing output for each one.
        for (var hammerEntry : hammerMap.entrySet()) {
            String fid = hammerEntry.getValue().getFid();
            Feature feat = repGenome.getFeature(fid);
            String roleId;
            if (feat == null)
                roleId = "(null)";
            else {
                var roles = this.getFeatureRoles(feat);
                if (roles.size() < 1)
                    roleId = "";
                else
                    roleId = roles.stream().map(x -> x.getId()).collect(Collectors.joining(", "));
            }
            this.printFields(hammerEntry.getKey(), fid, repId, repName, roleId, feat.getPegFunction());
        }
    }

    @Override
    protected void finishReport() {
    }

}
