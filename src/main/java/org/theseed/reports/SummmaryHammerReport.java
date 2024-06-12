/**
 *
 */
package org.theseed.reports;

import java.io.IOException;
import java.util.Map;

import org.theseed.basic.ParseFailureException;
import org.theseed.proteins.hammer.HammerDb.Source;

/**
 * This report outputs a count of the number of hammers per genome.
 *
 * @author Bruce Parrello
 *
 */
public class SummmaryHammerReport extends HammerReport {

    /**
     * Construct a summary hammer report.
     *
     * @param processor		controlling command processor
     *
     * @throws ParseFailureException
     * @throws IOException
     */
    public SummmaryHammerReport(IParms processor) throws ParseFailureException, IOException {
        super(processor);
    }

    @Override
    protected void initReport() {
        this.printFields("rep_id", "rep_name", "count");
    }

    @Override
    public void processGenome(String repId, String repName, Map<String, Source> hammerMap) {
        int count = hammerMap.size();
        this.printFields(repId, repName, Integer.toString(count));
    }

    @Override
    protected void finishReport() {
    }

}
