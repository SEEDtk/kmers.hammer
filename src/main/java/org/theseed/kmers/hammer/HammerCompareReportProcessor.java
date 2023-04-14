/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;

import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.hammer.HashHammerDb;
import org.theseed.utils.BaseReportProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command compares two hammer databases and outputs hammers found only in one or the other.  The positional
 * parameters should be the two hammer load files.  The report will be produced on the standard output.  If the
 * hammer databases are not extremely similar, the output can grow quite large.
 *
 * Note that hammer databases, because of their enormous size, cannot be handled by the standard join and compare
 * utilities.  (This is also why the output can grow quite large.)
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output file for report (if not STDOUT)
 *
 * @author Bruce Parrello
 *
 */
public class HammerCompareReportProcessor extends BaseReportProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(HammerCompareReportProcessor.class);

    // COMMAND-LINE OPTIONS

    /** file name of first hammer database */
    @Argument(index = 0, metaVar = "hammers1.tbl", usage = "load file for first hammer DB", required = true)
    private File hammer1File;

    /** file name of second hammer database */
    @Argument(index = 1, metaVar = "hammers2.tbl", usage = "load file for second hammer DB", required = true)
    private File hammer2File;

    @Override
    protected void setReporterDefaults() {
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        // Insure both hammer databases are readable.
        if (! hammer1File.canRead())
            throw new FileNotFoundException("First hammer file " + hammer1File + " is not found or unreadable.");
        if (! hammer2File.canRead())
            throw new FileNotFoundException("Second hammer file " + hammer2File + " is not found or unreadable.");
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // Write the report header.
        writer.println("file\thammer\tfid\tstrength");
        // We load one file into memory and then check the other file.
        long miss1 = this.compare(hammer1File, hammer2File, writer);
        long miss2 = this.compare(hammer2File, hammer1File, writer);
        log.info("{} hammers missing from {}, {} missing from {}.", miss1, hammer1File, miss2, hammer2File);
    }

    /**
     * Load the first hammer file into memory and check the second hammer file against it.
     *
     * @param hammers1		file of hammers to load into memory
     * @param hammers2		file of hammers to check
     * @param writer		output report writer
     *
     * @return the number of hammers in the second file not found in the first
     *
     * @throws ParseFailureException
     * @throws IOException
     */
    private long compare(File hammers1, File hammers2, PrintWriter writer) throws IOException, ParseFailureException {
        // This will count the missing hammers.
        long retVal = 0;
        // This will be displayed as the file name.
        String name = hammers2.getName();
        // Load the first file.
        log.info("Loading hammers from {}.", hammers1);
        HashHammerDb hammers = new HashHammerDb(hammers1);
        // Scan the second file.
        log.info("Checking hammers from {}.", hammers2);
        try (TabbedLineReader hammerStream = new TabbedLineReader(hammers2)) {
            // This counts the hammers checked.
            long count = 0;
            // This tracks the time between trace messages.
            long start = System.currentTimeMillis();
            // Loop through the input hammers.
            for (var line : hammerStream) {
                String hammer = line.get(0);
                count++;
                if (hammers.getSource(hammer) == null) {
                    // Here the hammer was not found.
                    retVal++;
                    writer.write(name + "\t" + hammer + "\t" + line.get(1) + "\t" + Double.toString(line.getDouble(2)));
                }
                if (log.isInfoEnabled() && System.currentTimeMillis() - start > 5000) {
                    start = System.currentTimeMillis();
                    log.info("{} hammers read from {}.  {} not found.", count, hammers2, retVal);
                }
            }
        }
        return retVal;
    }

}
