/**
 *
 */
package org.theseed.proteins.hammer;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.time.Duration;
import java.util.Collection;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.counters.CountMap;
import org.theseed.genome.Feature;
import org.theseed.io.TabbedLineReader;
import org.theseed.java.erdb.DbConnection;
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
         * @param genome	target genome
         * @param hammer	hammer that identifies it
         */
        public void updateHammerMap(String genome, String hammer) throws SQLException, IOException;

        /**
         * Create an empty hammer map.
         *
         * @param inFile	input file containing the hammers, for size estimates
         */
        public void createEmptyMap(File inFile) throws SQLException, IOException;

        public void close();

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
                    String genome = Feature.genomeOf(line.get(1));
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
                    loader.updateHammerMap(genome, hammer);
                    hCount++;
                    // The loads operate at very different speeds.  We update the log once per second.
                    if (log.isInfoEnabled() && System.currentTimeMillis() - logPoint >= 1000) {
                        log.info("{} hammers read.", hCount);
                        logPoint = System.currentTimeMillis();
                    }
                }
                if (this.kmerSize == 0)
                    throw new ParseFailureException("No hammers found in file " + inFile + ".");
            } catch (SQLException e) {
                // Convert to an IO exception for compatability.
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
        return this.findClosestInternal(seqs, this.kmerSize);
    }

    /**
     * Compute the closest genomes to a set of sequences.  This is the internal version of the
     * main method that facilitates access to the kmer size by the subclass.
     *
     * @param seqs		list of sequences to scan
     * @param kSize		kmer size to use
     *
     * @return a count map detailing the number of hits for each genome
     */
    protected abstract CountMap<String> findClosestInternal(Collection<Sequence> seqs, int kSize);

}
