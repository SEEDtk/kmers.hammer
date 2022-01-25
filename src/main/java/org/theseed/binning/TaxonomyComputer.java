/**
 *
 */
package org.theseed.binning;

import java.io.File;
import java.io.IOException;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This is the base class for taxonomy computation methods.  The basic task is to
 * take in a contigs FASTA file and spit out a taxonomic ID and name.
 *
 * @author Bruce Parrello
 *
 */
public abstract class TaxonomyComputer {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(TaxonomyComputer.class);
    /** dummy to return when the taxonomic class cannot be computed */
    public static final Result UNKNOWN_TAXON = new Result("6666666", "Unknown Prokaryote");
    /** dummy to use for empty bins */
    public static final Result VIRTUAL_GENOME = new Result("", "");

    /**
     * This object is used to pass back the results of the computation.
     */
    public static class Result {

        /** taxonomic ID */
        private String id;
        /** scientific name */
        private String name;

        /**
         * Construct a new taxonomic result object.
         *
         * @param taxon		taxonomic id
         * @param sci_name	scientific name
         */
        public Result(String taxon, String sci_name) {
            this.id = taxon;
            this.name = sci_name;
        }

        /**
         * @return the taxonomic ID
         */
        public String getId() {
            return this.id;
        }

        /**
         * @return the scientific name
         */
        public String getName() {
            return this.name;
        }

        @Override
        public int hashCode() {
            final int prime = 31;
            int result = 1;
            result = prime * result + ((this.id == null) ? 0 : this.id.hashCode());
            return result;
        }

        @Override
        public boolean equals(Object obj) {
            if (this == obj) {
                return true;
            }
            if (!(obj instanceof Result)) {
                return false;
            }
            Result other = (Result) obj;
            if (this.id == null) {
                if (other.id != null) {
                    return false;
                }
            } else if (!this.id.equals(other.id)) {
                return false;
            }
            return true;
        }

    }

    /**
     * Estimate the taxonomic classification of a genome stored in a contigs FASTA file.
     *
     * @param fastaFile		FASTA file to analyze
     *
     * @return the taxonomic classification found
     *
     * @throws IOException
     */
    public abstract Result analyzeFasta(File fastaFile) throws IOException;

}
