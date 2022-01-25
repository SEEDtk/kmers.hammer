/**
 *
 */
package org.theseed.binning;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.sequence.DnaInputStream;
import org.theseed.sequence.blast.BlastDB;
import org.theseed.sequence.blast.BlastHit;
import org.theseed.sequence.blast.BlastParms;
import org.theseed.sequence.blast.DnaBlastDB;

/**
 * This taxonomy computer accepts as input a FASTA BLAST database of seed proteins.  The comment should
 * contain a genome ID, a tab, and a scientific name.  The taxonomic ID is taken from the first part of
 * the genome ID, and the name from the scientific name.
 *
 * @author Bruce Parrello
 *
 */
public class BlastTaxonomyComputer extends TaxonomyComputer {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(BlastTaxonomyComputer.class);
    /** seed protein blast database */
    private BlastDB seedDb;
    /** BLAST parameters to use */
    private BlastParms parms;
    /** match pattern for extracting taxon ID and name out of a subject definition */
    private static final Pattern COMMENT_PATTERN = Pattern.compile("(\\d+)\\.\\d+\\s+(.+)");

    /**
     * Construct a BLAST-based taxonomy computer.
     */
    public BlastTaxonomyComputer(File seedFastaFile, BlastParms parms) throws IOException {
        // Note we only want the best match for each contig.
        this.parms = parms.clone().maxPerQuery(1);
        try {
            this.seedDb = DnaBlastDB.createOrLoad(seedFastaFile, 11);
        } catch (InterruptedException e) {
            throw new IOException("Error creating BLAST database: " + e.toString());
        }
    }

    @Override
    public Result analyzeFasta(File fastaFile) throws IOException {
        // Get the BLAST results.
        var fileStream = new DnaInputStream(fastaFile, 11);
        List<BlastHit> results = this.seedDb.blast(fileStream, this.parms);
        Result retVal;
        switch (results.size()) {
        case 0:
            retVal = TaxonomyComputer.UNKNOWN_TAXON;
            break;
        case 1:
            retVal = BlastTaxonomyComputer.getResult(results.get(0));
            break;
        default:
            // Here we must pick the best result.
            BlastHit best = results.get(0);
            int bestScore = best.getNumIdentical();
            for (int i = 1; i < results.size(); i++) {
                BlastHit result = results.get(i);
                int score = result.getNumIdentical();
                if (score > bestScore) {
                    bestScore = score;
                    best = result;
                }
            }
            retVal = BlastTaxonomyComputer.getResult(best);
        }
        return retVal;
    }

    /**
     * @return a taxonomic result object corresponding to the specified BLAST hit
     *
     * @param best		blast hit to convert
     *
     * @throws IOException
     */
    private static Result getResult(BlastHit best) throws IOException {
        Result retVal;
        String comment = best.getSubjectDef();
        Matcher m = COMMENT_PATTERN.matcher(comment);
        if (! m.matches())
            throw new IOException("Invalid seed protein FASTA file.  Comment for " + best.getSubjectId()
                    + " is badly formatted.");
        else
            retVal = new Result(m.group(1), m.group(2));
        return retVal;
    }

}
