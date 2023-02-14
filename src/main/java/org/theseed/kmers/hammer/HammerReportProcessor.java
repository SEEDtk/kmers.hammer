/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Map;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.hammer.HammerDb;
import org.theseed.reports.HammerReport;
import org.theseed.utils.BaseHammerUsageProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command produces a report on the hammers.  The positional parameter is the name of the representative
 * genome stats file for the repgen set from which the hammers were generated.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output file (if not STDOUT)
 * -b	batch size for queries
 *
 * --format			format of the report (default GENOMES)
 * --hType			type of hammer database (default MEMORY)
 * --method			voting method to use (default COUNT)
 * --file			file containing hammer database (either SQLite database or hammer flat file)
 * --url			URL of database (host and name, MySQL only)
 * --parms			database connection parameter string (MySQL only)
 * --type			database engine type
 * --roles			role definition table for computing roles
 * --reps			genome source containing representative genomes used to generate the hammers
 * --rSourceType	type of genome source for representative genomes (default DIR)
 *
 * @author Bruce Parrello
 *
 */
public class HammerReportProcessor extends BaseHammerUsageProcessor implements HammerReport.IParms {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(HammerReportProcessor.class);
    /** report writer */
    private HammerReport reportWriter;
    /** map of representative genome IDs to names */
    private Map<String, String> repGenomeMap;
    /** representative-genome source (if any) */
    private GenomeSource repGenomes;

    // COMMAND-LINE OPTIONS

    /** output report format */
    @Option(name = "--format", usage = "output report format")
    private HammerReport.Type reportType;

    /** representative genome source */
    @Option(name = "--reps", metaVar = "P3Eval/GTOXXX", usage = "directory or file containing representative genomes")
    private File repDir;

    /** type of representative genome source */
    @Option(name = "--rSourceType", usage = "type of representative genome source")
    private GenomeSource.Type rSourceType;

    /** role definition file */
    @Option(name = "--roles", metaVar = "roles.in.subsystems", usage = "role definition file")
    private File roleFile;

    /** representative-genome stats file for this hammer set */
    @Argument(index = 0, metaVar = "repXXX.stats.tbl", usage = "representative-genome statistics file from which the hammers were built")
    private File repStatsFile;

    @Override
    protected void setHammerDefaults() {
        this.reportType = HammerReport.Type.GENOMES;
        this.repDir = null;
        this.rSourceType = GenomeSource.Type.DIR;
        this.roleFile = null;
    }

    @Override
    protected void validateHammerParms() throws IOException, ParseFailureException {
        // Load the representative genome IDs and names.
        if (! this.repStatsFile.canRead())
            throw new FileNotFoundException("Repgen stats file " + this.repStatsFile + " is not found or unreadable.");
        this.repGenomeMap = TabbedLineReader.readMap(this.repStatsFile, "rep_id", "rep_name");
        log.info("{} representative genomes found in {}.", this.repGenomeMap.size(), this.repStatsFile);
        // If there is a representative-genome source, set it up.
        if (this.repDir == null) {
            this.repGenomes = null;
            log.info("No representative-genome source provided.");
        } else {
            this.repGenomes = this.rSourceType.create(this.repDir);
            log.info("{} representative genomes found in {}.", this.repGenomes.size(), this.repDir);
        }
        // Create the report.  This also validates the report's tuning parameters.
        log.info("Initializing report of type {}.", this.reportType);
        this.reportWriter = this.reportType.create(this);
    }

    @Override
    protected void runHammers(HammerDb hammerDb, PrintWriter writer) throws Exception {
        // Initialize the report.
        this.reportWriter.openReport(writer);
        // Loop through the representative genomes.
        int repCount = 0;
        for (var repEntry : this.repGenomeMap.entrySet()) {
            repCount++;
            String repId = repEntry.getKey();
            log.info("Processing representative {} of {}: {} ({})", repCount, this.repGenomeMap.size(), repId, repEntry.getValue());
            Map<String, HammerDb.Source> hammerMap = hammerDb.findGenomeHammers(repId);
            log.info("{} hammers found for {}.", hammerMap.size(), repId);
            this.reportWriter.processGenome(repId, repEntry.getValue(), hammerMap);
            writer.flush();
        }
        // Finish the report.
        this.reportWriter.closeReport();
    }

    @Override
    public File getRoleFile() {
        return this.roleFile;
    }

    @Override
    public GenomeSource getRepGenomes() {
        return this.repGenomes;
    }

}
