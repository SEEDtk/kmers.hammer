/**
 *
 */
package org.theseed.reports;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;
import org.theseed.proteins.hammer.HammerDb.Source;
import org.theseed.utils.ParseFailureException;

/**
 * This is the base class for all hammer reports.
 *
 * @author Bruce Parrello
 *
 */
public abstract class HammerReport {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(HammerReport.class);
    /** output print writer */
    private PrintWriter writer;
    /** representative-genome source (if any) */
    private GenomeSource repGenomes;
    /** role definition map (if any) */
    private RoleMap roles;

    /**
     * This interface defines the methods that must be supported for a controlling command processor.
     * It is used by the report writers to access the necessary tuning parameters.
     */
    public interface IParms {

        /**
         * @return the name of the role definition file, or NULL if there is none
         */
        public File getRoleFile();

        /**
         * @return the genome source for the representative genomes, or NULL if there is none
         */
        public GenomeSource getRepGenomes();

    }

    /**
     * This enumeration defines the different report types.
     */
    public static enum Type {
        /** hammers organized by genome, showing the roles */
        GENOMES {
            @Override
            public HammerReport create(IParms processor) throws ParseFailureException, IOException {
                return new GenomeHammerReport(processor);
            }
        };

        /**
         * @return a hammer report writer of this type
         *
         * @param processor		controlling command processor
         */
        public abstract HammerReport create(IParms processor) throws ParseFailureException, IOException;

    }

    /**
     * Construct a hammer report for a specified command processor.
     *
     * @param processor		controlling command processor
     *
     * @throws ParseFailureException
     * @throws IOException
     */
    public HammerReport(IParms processor) throws ParseFailureException, IOException {
        // Get the role definitions (if any).
        File roleFile = processor.getRoleFile();
        if (roleFile == null)
            this.roles = null;
        else {
            log.info("Loading role definitions from {}.", roleFile);
            this.roles = RoleMap.load(roleFile);
        }
        // Get the representative-genome source (if any).
        this.repGenomes = processor.getRepGenomes();
    }

    /**
     * Initialize the report.
     *
     * @param repWriter		output stream for the report
     */
    public void openReport(PrintWriter repWriter) {
        this.writer = repWriter;
        this.initReport();

    }

    /**
     * Initialize the report for the subclass.  Usually, this means writing the header line and setting
     * up data structures.
     */
    protected abstract void initReport();


    /**
     * Write a line of output.
     *
     * @param line		output line to writer
     */
    public void println(String line) {
        this.writer.println(line);
    }

    /**
     * Write a set of output fields.
     *
     * @param strings	array of strings containing the field values
     */
    public void printFields(String... strings) {
        this.writer.println(StringUtils.join(strings, '\t'));
    }

    /**
     * Process the hammers for a single representative genome.
     *
     * @param repId			representative genome ID
     * @param repName		representative genome name
     * @param hammerMap		map of hammers to source descriptors for each hammer generated from the genome
     */
    public abstract void processGenome(String repId, String value, Map<String, Source> hammerMap);

    /**
     * Summarize and complete the report.
     */
    public void closeReport() {
        this.finishReport();
    }

    /**
     * Finish the report.  Reports that contain only summary information are written in their entirety here.
     */
    protected abstract void finishReport();

    /**
     * @return the identified representative genome
     *
     * @param repId		ID of the desired genome
     */
    public Genome getRepGenome(String repId) {
        return this.repGenomes.getGenome(repId);
    }

    /**
     * @return the roles for a genome feature
     *
     * @param feat	feature whose role descriptors are desired
     */
    public List<Role> getFeatureRoles(Feature feat) {
        var retVal = feat.getUsefulRoles(this.roles);
        return retVal;
    }

    /**
     * @return TRUE if a representative-genome source is present
     */
    public boolean hasRepGenomes() {
        return this.repGenomes != null;
    }

    /**
     * @return TRUE if a role definition map is present
     */
    public boolean hasRoles() {
        return this.roles != null;
    }

}
