/**
 *
 */
package org.theseed.reports;

import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.theseed.basic.ParseFailureException;
import org.theseed.counters.WeightMap;
import org.theseed.genome.Genome;
import org.theseed.proteins.hammer.HammerDb;
import org.theseed.proteins.hammer.HammerDb.Source;
import org.theseed.sequence.Sequence;

/**
 * This report scans all the hammers in a representative genome and insures they belong to its subset of the
 * hammer database.
 *
 * @author Bruce Parrello
 *
 */
public class ValidationHammerReport extends HammerReport {

    /**
     * Create a new hammer validation report.
     *
     * @param processor		controlling command processor
     *
     * @throws ParseFailureException
     * @throws IOException
     */
    public ValidationHammerReport(IParms processor) throws ParseFailureException, IOException {
        super(processor);
        if (! this.hasRepGenomes())
            throw new ParseFailureException("Representative-genome source is required for this type of report.");
    }

    @Override
    protected void initReport() {
        // Write the report header.
        this.printFields("rep_id", "rep_name", "base pairs", "expected", "found", "alien", "total", "best", "best_score",
                "num_aliens", "max_alien", "mean_alien");
    }

    @Override
    public void processGenome(String repId, String value, Map<String, Source> hammerMap) {
        // Process each contig against the main hammer map.  We need to know how many of the
        // expected hammers are found and how many alien hammers were found.
        HammerDb hammers = this.getHammers();
        // Get the representative genome.
        Genome genome = this.getRepGenome(repId);
        log.info("Processing genome {} with {} expected hammers.", genome, hammerMap.size());
        // Get all the contig sequences.
        Collection<Sequence> contigs = genome.getContigs().stream().map(x -> new Sequence(x.getId(), "", x.getSequence()))
                .collect(Collectors.toList());
        // Loop through the contigs, accumulating the hammers found.
        Set<String> hammerSet = contigs.stream().flatMap(x -> hammers.findHammers(x.getSequence()).stream())
                .collect(Collectors.toSet());
        log.info("{} hammers found in genome.", hammerSet.size());
        // Count the number of expected hammers found in the genome and compute the number of unexpected ones.
        int found = (int) hammerMap.keySet().stream().filter(x -> hammerSet.contains(x)).count();
        int alien = hammerSet.size() - found;
        // Now score the genome.
        WeightMap scoreMap = hammers.findClosest(contigs);
        // Show the top-scoring genome.
        WeightMap.Count best = scoreMap.getBestEntry();
        // Compute the mean and the max of the others.
        double maxWeight = 0.0;
        double totalWeight = 0.0;
        int numWeights = 0;
        for (var count : scoreMap.counts()) {
            if (! repId.contentEquals(count.getKey())) {
                double weight = count.getCount();
                if (weight > maxWeight)
                    maxWeight = weight;
                totalWeight += weight;
                numWeights++;
            }
        }
        String mean = "";
        if (numWeights > 0)
            mean = Double.toString(totalWeight / numWeights);
        String max = "";
        if (maxWeight > 0.0)
            max = Double.toString(maxWeight);
        this.printFields(genome.getId(), genome.getName(), Integer.toString(genome.getLength()),
                Integer.toString(hammerMap.size()), Integer.toString(found), Integer.toString(alien),
                Integer.toString(hammerSet.size()), best.getKey(), Double.toString(best.getCount()),
                Integer.toString(numWeights), max, mean);
    }

    @Override
    protected void finishReport() {
    }

}
