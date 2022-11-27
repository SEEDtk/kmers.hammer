/**
 *
 */
package org.theseed.proteins.hammer;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.time.Duration;
import java.util.Collection;
import java.util.List;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.counters.CountMap;
import org.theseed.genome.Feature;
import org.theseed.io.TabbedLineReader;
import org.theseed.java.erdb.DbConnection;
import org.theseed.locations.Location;
import org.theseed.sequence.Sequence;
import org.theseed.utils.ParseFailureException;

/**
 * This object represents a thor-hammer database.  Such databases are huge and take up a great deal of
 * memory; however, the code to manage them is small and light.
 *
 * A thor-hammer database maps DNA kmers of a specified length to representative-genome IDs.  If
 * given a list of sequences, it will compute the closest representative genomes to that sequence
 * list.
 *
 * This is the base class.  There are several subclasses that store the hammers in different ways.
 *
 * @author Bruce Parrello
 *
 */
public abstract class HammerDb {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(HashHammerDb.class);
    /** kmer size to use */
    private int kmerSize;

    /**
     * This describes the interface that command processors must support to use the
     * HammerDb subclasses.
     */
    public interface IParms {

        /**
         * @return the database file
         */
        public File getDbFile();

        /**
         * @return the database connection for the hammer database
         */
        public DbConnection getConnection();

        /**
         * @return the batch size to use for database queries
         */
        public int getBatchSize();

    }

    /**
     * This describes the interface that a loading object must support to load the hammer
     * database.
     */
    public interface ILoader extends AutoCloseable {

        /**
         * Forge a connection between a hammer and a genome.
         *
         * @param fid		feature ID of the universal protein in the target genome
         * @param hammer	hammer that identifies it
         */
        public void updateHammerMap(String fid, String hammer) throws SQLException, IOException;

        /**
         * Create an empty hammer map.
         *
         * @param inFile	input file containing the hammers, for size estimates
         */
        public void createEmptyMap(File inFile) throws SQLException, IOException;

        public void close();

    }

    /**
     * This object describes a hammer hit.  Hammer hits are sorted by hit location and then feature ID, so they
     * are organized according to the contig hit.
     */
    public static class Hit implements Comparable<Hit> {

        /** contig location hit */
        private Location loc;
        /** feature ID of the hammer hit */
        private String fid;

        /**
         * Construct a hit descriptor.
         *
         * @param contig	contig ID
         * @param len		contig len
         * @param idx		0-based index of the kmer on the contig
         * @param dir		TRUE if the hit was on the plus strand, else FALSE
         * @param fid		ID of the feature from which the hammer was harvested
         * @param kSize		kmer size of the hammer
         *
         */
        protected Hit(String contig, int len, int idx, boolean dir, String fid, int kSize) {
            // We will save the start and end locations of the hit in here.
            int start;
            int end;
            // Initialize them according to the direction.
            if (dir) {
                // Here we have a hit on the + strand.  The locations are 1-based, so we adjust the 0-based idx.
                start = idx + 1;
                end = idx + kSize - 2;
            } else {
                // Here we have a hit on the - strand.  The start location is therefore counted from the end of the contig.
                start = len - idx;
                end = start - kSize + 1;
            }
            this.loc = Location.create(contig, start, end);
            this.fid = fid;
        }

        /**
         * @return the location of the hit
         */
        public Location getLoc() {
            return this.loc;
        }

        /**
         * @return the ID of the feature from which the hammer was harvested
         */
        public String getFid() {
            return this.fid;
        }

        /**
         * Specify the ID of the feature from which the hammer was harvested.
         *
         * @param fid	proposed new feature ID
         */
        public void setFid(String fid) {
            this.fid = fid;
        }

        /**
         * @return the ID of the genome from which the hammer was harvested
         */
        public String getGenomeId() {
            return Feature.genomeOf(this.fid);
        }

        @Override
        public int hashCode() {
            final int prime = 31;
            int result = 1;
            result = prime * result + ((this.fid == null) ? 0 : this.fid.hashCode());
            result = prime * result + ((this.loc == null) ? 0 : this.loc.hashCode());
            return result;
        }

        @Override
        public boolean equals(Object obj) {
            if (this == obj) {
                return true;
            }
            if (!(obj instanceof Hit)) {
                return false;
            }
            Hit other = (Hit) obj;
            if (this.fid == null) {
                if (other.fid != null) {
                    return false;
                }
            } else if (!this.fid.equals(other.fid)) {
                return false;
            }
            if (this.loc == null) {
                if (other.loc != null) {
                    return false;
                }
            } else if (!this.loc.equals(other.loc)) {
                return false;
            }
            return true;
        }

        @Override
        public int compareTo(Hit o) {
            int retVal = this.loc.compareTo(o.loc);
            if (retVal == 0)
                retVal = this.fid.compareTo(o.fid);
            return retVal;
        }

    }


    /**
     * This enumeration describes the types of hammer databases.
     */
    public static enum Type {
        MEMORY {
            @Override
            public HammerDb create(IParms processor) throws ParseFailureException, IOException {
                return new HashHammerDb(processor);
            }
        }, ERDB {
            @Override
            public HammerDb create(IParms processor) {
                return new SqlHammerDb(processor);
            }
        };

        /**
         * @return a hammer database of this type
         *
         * @param processor		controlling command processor
         *
         * @throws IOException
         * @throws ParseFailureException
         * @throws SQLException
         */
        public abstract HammerDb create(IParms processor) throws IOException, ParseFailureException, SQLException;
    }

    /**
     * Load a hammer database from a hammer input file.
     *
     * @param inFile		input file containing the hammers
     *
     * @throws ParseFailureException
     * @throws IOException
     */
    public void load(File inFile) throws ParseFailureException, IOException {
        try (TabbedLineReader inStream = new TabbedLineReader(inFile)) {
            // Denote we do not yet know the kmer size.
            this.kmerSize = 0;
            // Start the timers and the counting.
            int hCount = 0;
            long start = System.currentTimeMillis();
            long logPoint = start;
            // Get the loader object.
            try (ILoader loader = this.getLoader()) {
                // Create the map.
                loader.createEmptyMap(inFile);
                // Read in the hammers.
                for (TabbedLineReader.Line line : inStream) {
                    String fid = line.get(1);
                    String hammer = line.get(0);
                    if (this.kmerSize != hammer.length()) {
                        if (this.kmerSize == 0)
                            this.kmerSize = hammer.length();
                        else {
                            log.error("Invalid kmer \"{}\" in hammer file {}: length is not {}.",
                                    hammer, inFile, this.kmerSize);
                            throw new ParseFailureException("Invalid kmer \"" + hammer + "\" in hammer file (bad length).");
                        }
                    }
                    // Add this hammer to the map.
                    loader.updateHammerMap(fid, hammer);
                    hCount++;
                    // The loads operate at very different speeds.  We update the log once per 5 seconds.
                    if (log.isInfoEnabled() && System.currentTimeMillis() - logPoint >= 5000) {
                        log.info("{} hammers read.", hCount);
                        logPoint = System.currentTimeMillis();
                    }
                }
                if (this.kmerSize == 0)
                    throw new ParseFailureException("No hammers found in file " + inFile + ".");
            } catch (SQLException e) {
                // Convert to an IO exception for compatibility.
                throw new IOException("SQL Error: " + e.toString());
            }
            if (log.isInfoEnabled()) {
                Duration length = Duration.ofMillis(System.currentTimeMillis() - start);
                log.info("Hammer map contains {} hammers. {} to load.", hCount, length.toString());
            }
        }
    }

    /**
     * The loading object contains resources that have to be freed at the end of loading.
     * It is also
     * @return a loading object
     */
    protected abstract ILoader getLoader();

    /**
     * Set the kmer size for the database.
     *
     * @param newSize	new size to set
     */
    protected void setKmerSize(int newSize) {
        this.kmerSize = newSize;
    }

    /**
     * Compute the closest genomes to a set of sequences.  We accumulate the hammer hits in a count
     * map keyed by genome ID and simply return the count map to the caller.
     *
     * This method must be thread-safe, so it does not alter any object fields.
     *
     * @param seqs	list of sequences to scan
     *
     * @return a count map detailing the number of hits for each genome
     */
    public CountMap<String> findClosest(Collection<Sequence> seqs) {
        var retVal = new CountMap<String>();
        // Get the reverse complements.
        List<Sequence> seqList = reverseAll(seqs);
        // Add the original sequences.
        seqList.addAll(seqs);
        // Find the closest genomes to the set of sequences.
        this.findClosestInternal(retVal, seqList, this.kmerSize);
        // Return the results.
        return retVal;
    }

    /**
     * Reverse all the sequences in a sequence collection.
     *
     * @param seqs		collection of sequences to process
     *
     * @return a collection of new sequences containing the reverse complements
     */
    public static List<Sequence> reverseAll(Collection<Sequence> seqs) {
        return seqs.stream().map(x -> x.reverse()).collect(Collectors.toList());
    }

    /**
     * @return the individual hits for a collection of sequences
     */
    public Collection<Hit> findHits(Collection<Sequence> seqs) {
        // Process the forward direction.
        TreeSet<Hit> retVal = new TreeSet<Hit>();
        this.findHitsInternal(retVal, seqs, this.kmerSize, true);
        // Get the reverse complements.
        List<Sequence> revs = reverseAll(seqs);
        this.findHitsInternal(retVal, revs, this.kmerSize, false);
        // Return the full set of hits.
        return retVal;
    }

    /**
     * Extract the hammer hits from a collection of sequences.  There is no need to check reverse complements:  this
     * is facilitated by the framework.
     *
     * @param collection	output collection for the hits
     * @param seqs			collection of sequences to process
     * @param kSize			kmer size
     * @param dir			TRUE if this is a collection of plus strands, else FALSE
     */
    protected abstract void findHitsInternal(Collection<Hit> collection, Collection<Sequence> seqs, int kSize, boolean dir);

    /**
     * Compute the closest genomes to a set of sequences.  This is the internal version of the
     * main method that facilitates access to the kmer size by the subclass.  There is no need to check reverse
     * complements:  this is facilitated by the framework.
     *
     * @param map		count map to contain the results
     * @param seqs		list of sequences to scan
     * @param kSize		kmer size to use
     *
     * @return a count map detailing the number of hits for each genome
     */
    protected abstract void findClosestInternal(CountMap<String> map, Collection<Sequence> seqs, int kSize);

}
