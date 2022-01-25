/**
 *
 */
package org.theseed.binning;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.Sequence;

/**
 * This object manages the creation of bin FASTA files.  It takes as parameters an output directory and
 * a method of deciding how contigs are sorted into bins.  The client then passes it contig sequences that
 * are sorted into separate bin files created by this object.
 *
 * NOTE that the bin ID must be valid as a part of a file name.
 *
 * @author Bruce Parrello
 *
 */
public class BinBuilder implements AutoCloseable {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(BinBuilder.class);
    /** binning rule */
    private IRule rule;
    /** map of bin IDs to control objects */
    private Map<String, Stats> binMap;
    /** bin control object for rejected contigs */
    private Stats rejected;
    /** default estimated coverage */
    private static final double FALLBACK_COVERAGE = 50.0;
    /** coverage match pattern for comment */
    private static final Pattern COMMENT_COVERAGE_PATTERN = Pattern.compile("\\b(?:covg|coverage|cov|multi)[= ]([0-9.]+)",
            Pattern.CASE_INSENSITIVE);
    /** coverage match pattern for label */
    private static final Pattern LABEL_COVERAGE_PATTERN = Pattern.compile("(?:coverage|covg|cov)_([0-9.]+)",
            Pattern.CASE_INSENSITIVE);

    /**
     * This interface must be supported by any valid binning rule object.
     */
    public interface IRule {

        /**
         * @return the bin IDs
         */
        public Set<String> getBinIds();

        /**
         * Determine the destination bin for a contig.
         *
         * @param seq	sequence of the contig (includes label, comment, and DNA)
         *
         * @return the ID of the target bin, or NULL if the contig should be rejected
         */
        public String choose(Sequence seq);

    }

    /**
     * This object describes a bin.
     */
    public static class Stats implements Comparable<Stats> {

        /** ID of the bin */
        private String id;
        /** number of contigs in the bin */
        private int num;
        /** sequence length of the bin */
        private long length;
        /** file name of the bin */
        private File location;
        /** output stream for the bin, or NULL if it is closed */
        private FastaOutputStream outStream;
        /** TRUE if this is the bin for unplaced contigs */
        private boolean virtual;

        /**
         * Construct a blank bin statistics object.
         *
         * @param binId		ID of the bin
         * @param outDir	output directory for the bin's FASTA file
         */
        protected Stats(String binId, File outDir) {
            this.id = binId;
            this.init();
            this.location = new File(outDir, "bin." + binId + ".fasta");
            this.virtual = false;
        }

        /**
         * Construct a bin statistics object for rejected contigs.
         */
        protected Stats(File outDir) {
            this.id = "(rejected)";
            this.init();
            this.location = new File(outDir, "unplaced.fasta");
            this.virtual = true;
        }

        /**
         * Initialize the counters and totals.
         */
        private void init() {
            this.num = 0;
            this.length = 0;
            this.outStream = null;
        }

        /**
         * Store a contig in this bin.
         *
         * @throws IOException
         */
        protected void store(Sequence seq) throws IOException {
            if (this.outStream == null) {
                // Here we must open the output stream so we can put the contig in it.
                log.info("Creating file {} for bin {}.", this.location, this.id);
                this.outStream = new FastaOutputStream(this.location);
            }
            // Put the contig in the FASTA file.
            this.outStream.write(seq);
            this.num++;
            this.length += seq.length();
        }

        /**
         * Close the stream for this bin.
         */
        protected void close() {
            if (this.outStream != null)
                this.outStream.close();
        }

        /**
         * @return the ID of the bin
         */
        public String getId() {
            return this.id;
        }

        /**
         * @return the number of contigs in the bin
         */
        public int getNum() {
            return this.num;
        }

        /**
         * @return the sequence length of the bin
         */
        public long getLength() {
            return this.length;
        }

        /**
         * @return the file location of the bin
         */
        public File getLocation() {
            return this.location;
        }

        @Override
        public int compareTo(Stats o) {
            // Sort virtual bins to the end.
            int retVal = Boolean.compare(o.virtual, this.virtual);
            if (retVal == 0) {
                // Put the biggest bins first.
                retVal = Long.compare(o.length, this.length);
                if (retVal == 0) {
                    // If all else fails, sort by ID.
                    retVal = this.id.compareTo(o.id);
                }
            }
            return retVal;
        }

        /**
         * @return TRUE if this is a virtual bin
         */
        public boolean isVirtual() {
            return this.virtual;
        }

    }

    /**
     * Construct a bin builder.
     *
     * @param rule		binning rule to use
     * @param dir		output directory
     */
    public BinBuilder(IRule rule, File dir) {
        this.rule = rule;
        // Get the list of bin IDs and build the bin map.
        var binIds = rule.getBinIds();
        this.binMap = binIds.stream().collect(Collectors.toMap(x -> x, x -> new Stats(x, dir)));
        // Set up the rejected-contigs bin.
        this.rejected = new Stats(dir);
    }

    /**
     * Bin a contig.
     *
     * @param seq	contig to bin
     *
     * @return the ID of the chosen bin, or NULL if the contig was rejected
     *
     * @throws IOException
     */
    public String binContig(Sequence seq) throws IOException {
        String retVal = this.rule.choose(seq);
        Stats bin = this.binMap.getOrDefault(retVal, this.rejected);
        bin.store(seq);
        return retVal;
    }

    /**
     * @return all the bin control objects for real bins
     */
    public Collection<Stats> getBins() {
        return this.binMap.values();
    }

    /**
     * @return all the nonempty bin control objects (including the virtual bin for rejected contigs)
     */
    public Collection<Stats> getNonEmptyBins() {
        List<Stats> retVal = new ArrayList<Stats>(this.binMap.size() + 1);
        // Get all the real bins.
        this.getBins().stream().filter(x -> x.length > 0).forEach(x -> retVal.add(x));
        // Add the fake bin for rejected contigs.
        if (this.rejected.length > 0)
            retVal.add(this.rejected);
        // Sort the results.
        Collections.sort(retVal);
        return retVal;
    }

    /**
     * @return the control object for rejected contigs
     */
    public Stats getRejected() {
        return this.rejected;
    }

    @Override
    public void close() throws Exception {
        // Close any open output streams.
        for (Stats stat : this.binMap.values()) {
            stat.close();
        }
    }

    /**
     * Compute the coverage for a contig.  We look for special keywords in the contig label and comment.
     *
     * @param seq	contig sequence object
     *
     * @return the estimated coverage
     */
    public static double computeCoverage(Sequence seq) {
        String coverage = null;
        // First, search the comment.
        String comment = seq.getComment();
        Matcher m = COMMENT_COVERAGE_PATTERN.matcher(comment);
        if (m.find())
            coverage = m.group(1);
        else {
            // Failing that, search the label.
            String label = seq.getLabel();
            m = LABEL_COVERAGE_PATTERN.matcher(label);
            if (m.find())
                coverage = m.group(1);
        }
        double retVal;
        if (coverage == null)
            retVal = FALLBACK_COVERAGE;
        else
            retVal = Double.valueOf(coverage);
        return retVal;
    }

    /**
     * Determine the contig's coverage and length and decide whether or not the contig is acceptable.
     *
     * @param seq		sequence object for the contig
     * @param minCovg	minimum acceptable coverage
     * @param minLen	minimum acceptable length
     *
     * @return TRUE if the contig is acceptable, else FALSE
     */
    public static boolean contigFilter(Sequence seq, double minCovg, int minLen) {
        boolean retVal = false;
        // Start with the length.
        int len = seq.length();
        if (len >= minLen) {
            // Now check the coverage.
            double covg = computeCoverage(seq);
            retVal = (covg >= minCovg);
        }
        return retVal;
    }

}
