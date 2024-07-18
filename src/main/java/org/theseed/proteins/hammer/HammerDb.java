/**
 *
 */
package org.theseed.proteins.hammer;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.time.Duration;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Contig;
import org.theseed.genome.Feature;
import org.theseed.io.TabbedLineReader;
import org.theseed.java.erdb.DbConnection;
import org.theseed.locations.Location;
import org.theseed.sequence.ISequence;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.fastq.SeqRead;

/**
 * This object represents a thor-hammer database.  Such databases are huge and take up a great deal of
 * memory; however, the code to manage them is small and light.
 *
 * A thor-hammer database maps DNA kmers of a specified length to representative-genome IDs.  If
 * given a list of sequences, it will compute the closest representative genomes to that sequence
 * list.
 *
 * Hammers are stored as long integers, which means the maximum hammer size is 31 characters.
 *
 * This is the base class.  There are several subclasses that store the hammers in different ways.
 *
 * @author Bruce Parrello
 *
 */
public abstract class HammerDb {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(HammerDb.class);
    /** kmer size to use */
    private int kmerSize;
    /** number of quality rejections */
    private int badQual;
    /** set of roles found */
    private Set<String> roles;
    /** method to use for counting hits */
    private Method countMethod = Method.COUNT;

    /**
     * This enumeration describes the counting method to be used for finding the closest match.
     */
    public static enum Method {
        STRENGTH {
            @Override
            public double getWeight(ISource hitSource) {
                return hitSource.getStrength();
            }
        }, COUNT {
            @Override
            public double getWeight(ISource hitSource) {
                return 1.0;
            }
        };

        /**
         * @return the weight to use when counting this hammer hit
         *
         * @param hitSource		hammer hit source
         */
        public abstract double getWeight(ISource hitSource);

    }

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
         * @param roleId	role ID of the universal protein
         * @param hammer	hammer that identifies the feature
         * @param str		strength of the hammer
         */
        public void updateHammerMap(String fid, String roleId, String hammer, double str) throws SQLException, IOException;

        /**
         * Create an empty hammer map.
         *
         * @param inFile	input file containing the hammers, for size estimates
         */
        public void createEmptyMap(File inFile) throws SQLException, IOException;

        public void close();

    }

    /**
     * This object describes a hammer's source information.  It includes the source feature ID, role ID, and
     * strength.
     */
    public static class Source implements ISource, HammerMap.IScore {

        /** source feature ID for the hammer */
        private String fid;
        /** source role ID for the hammer */
        private String roleId;
        /** strength of the hammer, a ranking from 0 to 1 */
        private double strength;

        /**
         * Construct a new hammer source descriptor.
         *
         * @param fid		source feature ID
         * @param roleId	role ID for the feature
         * @param str		strength of the hammer
         */
        public Source(String fid, String roleId, double str) {
            this.fid = fid;
            this.roleId = roleId;
            this.strength = str;
        }

        /**
         * @return the source feature ID
         */
        public String getFid() {
            return this.fid;
        }

        /**
         * @return the strength ratio
         */
        public double getStrength() {
            return this.strength;
        }

        /**
         * @return the ID of the hammer's source genome
         */
        public String getGenomeId() {
            return Feature.genomeOf(this.fid);
        }

        /**
         * @return the ID of the hammer's source role
         */
        public String getRole() {
            return this.roleId;
        }

        @Override
        public boolean isBadHammer() {
            return false;
        }

        @Override
        public void setBadHammer() {
        }

    }

    /**
     * This interface can be used for any object that contains a hammer's role ID and strength.
     */
    public interface ISource {

        /**
         * @return the ID of the hammer's source role
         */
        public String getRole();

        /**
         * @return the hammer's hit strength
         */
        public double getStrength();

    }

    /**
     * This object describes a hammer hit.  Hammer hits are sorted by hit location and then feature ID, so they
     * are organized according to the contig hit.
     */
    public static class Hit implements Comparable<Hit>, ISource {

        /** contig location hit */
        private Location loc;
        /** feature ID of the hammer hit */
        private String fid;
        /** role ID of the feature */
        private String roleId;
        /** strength of the hammer */
        private double strength;
        /** hammer string */
        private String hammer;

        /**
         * Construct a hit descriptor.
         *
         * @param hammer	hammer string
         * @param contig	contig ID
         * @param len		contig len
         * @param idx		0-based index of the kmer on the contig
         * @param dir		TRUE if the hit was on the plus strand, else FALSE
         * @param fid		ID of the feature from which the hammer was harvested
         * @param role		ID of the feature's role
         * @param kSize		kmer size of the hammer
         * @param str		strength of the hammer
         *
         */
        protected Hit(String hammer, String contig, int len, int idx, boolean dir, String fid, String role, int kSize, double str) {
            // We will save the start and end locations of the hit in here.
            int start;
            int end;
            // Initialize them according to the direction.
            if (dir) {
                // Here we have a hit on the + strand.  The locations are 1-based, so we adjust the 0-based idx.
                start = idx + 1;
                end = idx + kSize;
            } else {
                // Here we have a hit on the - strand.  The start location is therefore counted from the end of the contig.
                start = len - idx;
                end = start - kSize + 1;
            }
            this.loc = Location.create(contig, start, end);
            this.fid = fid;
            this.strength = str;
            this.roleId = role;
            this.hammer = hammer;
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

        /**
         * @return the strength of the hammer that made the hit
         */
        public double getStrength() {
            return this.strength;
        }

        /**
         * Update the source feature ID and strength.
         *
         * @param fid2		hammer source feature ID
         * @param role2		hammer source role ID
         * @param strength2	strength ratio
         */
        public void setHammerSource(String fid2, String role2, double strength2) {
            this.fid = fid2;
            this.roleId = role2;
            this.strength = strength2;
        }

        @Override
        public String getRole() {
            return this.roleId;
        }

        /**
         * @return the hammer string
         */
        public String getHammer() {
            return this.hammer;
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

            @Override
            public boolean isLoadFileKnown() {
                return true;
            }
        }, ERDB {
            @Override
            public HammerDb create(IParms processor) {
                return new SqlHammerDb(processor);
            }

            @Override
            public boolean isLoadFileKnown() {
                return false;
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

        /**
         * @return TRUE if this database type supports the get-load-file function
         */
        public abstract boolean isLoadFileKnown();
    }

    /**
     * Construct a blank, empty hammer database.
     */
    protected HammerDb() {
        this.badQual = 0;
        this.roles = new TreeSet<String>();
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
                // Read in the hammers.
                for (TabbedLineReader.Line line : inStream) {
                    String fid = line.get(1);
                    String hammer = line.get(0);
                    if (this.kmerSize != hammer.length()) {
                        if (this.kmerSize == 0) {
                            this.kmerSize = hammer.length();
                            loader.createEmptyMap(inFile);
                        } else {
                            log.error("Invalid kmer \"{}\" in hammer file {}: length is not {}.",
                                    hammer, inFile, this.kmerSize);
                            throw new ParseFailureException("Invalid kmer \"" + hammer + "\" in hammer file (bad length).");
                        }
                    }
                    double strength = line.getDouble(2);
                    String role = line.get(6);
                    // Add this hammer to the map.
                    loader.updateHammerMap(fid, role, hammer, strength);
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
     * Specify the scoring method for counting hits.
     *
     * @param method	new method to use
     */
    public void setMethod(Method method) {
        this.countMethod = method;
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
        if (newSize > 31)
            throw new IllegalArgumentException("Kmer size for hammers must be 31 or less.");
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
     * @return a weight map detailing the hit score for each genome
     */
    public ScoreMap findClosest(Collection<Sequence> seqs) {
        var retVal = new ScoreMap();
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
     * @return the individual hits for a collection of sequences, sorted by location
     */
    public SortedSet<Hit> findHits(Collection<Sequence> seqs) {
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
     * @return the set of hammers found in a sequence
     *
     * @param seq	sequence to search
     */
    public Set<String> findHammers(String seq) {
        var retVal = new HashSet<String>(100);
        this.findHammersInternal(retVal, seq, this.kmerSize);
        String rSeq = Contig.reverse(seq);
        this.findHammersInternal(retVal, rSeq, this.kmerSize);
        return retVal;
    }

    /**
     * Find all the hammers in the specified single-strand sequence and add them to the specified output set.
     *
     * @param hammerSet		output set for hammers found
     * @param seq			sequence to search
     * @param kSize			kmer size to use
     */
    protected abstract void findHammersInternal(HashSet<String> hammerSet, String seq, int kSize);

    /**
     * Find all the hammers for a specified source genome.
     *
     * @param genomeId		ID of the source genome
     *
     * @return a map from each hammer in the genome to its source
     */
    public abstract Map<String, Source> findGenomeHammers(String genomeId);

    /**
     * Count a hit in a result map.
     *
     * @param scoreMap		weight map containing the hit score for each genome
     * @param hitSource		hit to count
     * @param count			number of times to count the hit
     */
    protected void countHit(ScoreMap scoreMap, Source hitSource, int count) {
        scoreMap.count(hitSource.getGenomeId(), this.countMethod.getWeight(hitSource) * count, hitSource.getRole());
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
    protected abstract void findHitsInternal(Collection<Hit> collection, Collection<? extends ISequence> seqs, int kSize, boolean dir);

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
    protected abstract void findClosestInternal(ScoreMap map, Collection<Sequence> seqs, int kSize);

    /**
     * @return the name of the hammer load file, or NULL if it is not available
     */
    public abstract File getLoadFile();

    /**
     * Encode a hammer into a long integer.  The hammer must not have any ambiguity characters.
     *
     * @param hammer	hammer to encode
     *
     * @return a long integer representation of the hammer, or a negative value if the hammer is invalid
     */
    public long encode(String hammer) {
        return HammerMap.encode(hammer, this.kmerSize);
    }

    /**
     * Decode a hammer from a long integer.
     *
     * @param coded		long integer to decode
     *
     * @return the hammer sequence string
     */
    public String decode(long coded) {
        return HammerMap.decode(coded, this.kmerSize);
    }

    /**
     * @return the kmer size for this hammer database
     */
    public int getKmerSize() {
        return this.kmerSize;
    }

    /**
     * This is an optional method that gets the source for a hammer.  If it is not supported, it throws an
     * exception.
     *
     * @param hammer	hammer string to find
     *
     * @return the source for a hammer
     *
     * @throws ParseFailureException
     */
    public Source getSource(String hammer) throws ParseFailureException {
        throw new ParseFailureException("Get-source function not supported by this hammer database type.");
    }

    /**
     * @return the countMethod
     */
    public Method getCountMethod() {
        return this.countMethod;
    }

    /**
     * This method computes the hammer hits above a certain quality level.
     *
     * @param values	collection of sequences with quality data
     * @param minQual	minimum chance the hit is correct
     *
     * @return the set of hammer hits at or above the specified quality level
     */
    public SortedSet<Hit> findHits(Collection<SeqRead.Part> values, double minQual) {
        // We need a map of the parts and a list of the reversed parts.
        Map<String, SeqRead.Part> partMap = new HashMap<String, SeqRead.Part>(values.size() * 4 / 3 + 1);
        List<SeqRead.Part> revParts = new ArrayList<SeqRead.Part>(values.size());
        for (SeqRead.Part part : values) {
            partMap.put(part.getLabel(), part);
            SeqRead.Part rPart = SeqRead.Part.reverse(part);
            revParts.add(rPart);
        }
        // Get the hits for this collection of sequences in both directions.
        SortedSet<Hit> retVal = new TreeSet<Hit>();
        this.findHitsInternal(retVal, values, this.kmerSize, true);
        this.findHitsInternal(retVal, revParts, this.kmerSize, false);
        // Now loop through the hits, removing low-quality ones.
        Iterator<Hit> iter = retVal.iterator();
        while (iter.hasNext()) {
            Hit hit = iter.next();
            Location hitLoc = hit.getLoc();
            // Get the start of the hit location.  The length is always kmerSize.
            int hitPos = hitLoc.getLeft() - 1;
            // Get the sequence part hit.
            SeqRead.Part part = partMap.get(hitLoc.getContigId());
            // Compute the quality.
            String qString = part.getQual();
            double qual = SeqRead.qualChance(qString, hitPos, this.kmerSize);
            if (qual < minQual) {
                iter.remove();
                this.badQual++;
            }
        }
        return retVal;
    }

    /**
     * @return the number of bad-quality hits rejected
     */
    public int getBadQual() {
        return this.badQual;
    }

    /**
     * Save a role ID in the role set.
     *
     * @param roleId	role ID to save
     */
    protected void recordRole(String roleId) {
        this.roles.add(roleId);
    }

    /**
     * @return the set of roles in this hammer database
     */
    public Set<String> getRoles() {
        return this.roles;
    }

}
