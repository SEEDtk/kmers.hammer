/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.IOException;
import java.io.PrintWriter;

import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.io.TabbedLineReader;
import org.theseed.utils.BasePipeProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command take a hammer definition file in from the standard input and filters it by worthiness, producing
 * a smaller hammer definition file on the standard output.
 *
 * The positional parameter is the minimum worthiness level to keep.  It must be between 0 and 1.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more detailed log messages
 * -i	input file containing the original hammers (if not STDIN)
 * -o	output file to contain filtered hammers (if not STDOUT)
 *
 * @author Bruce Parrello
 *
 */
public class HammerWorthFilterProcessor extends BasePipeProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(HammerWorthFilterProcessor.class);
    /** worthiness column index */
    private int worthColIdx;

    // COMMAND-LINE OPTIONS

    /** minimum acceptable worthiness level */
    @Argument(index = 0, metaVar = "0.50", usage = "minimum worthiness level to keep", required = true)
    private double minWorth;

    @Override
    protected void setPipeDefaults() {
    }

    @Override
    protected void validatePipeParms() throws IOException, ParseFailureException {
        // Insure the worthiness is in range.
        if (this.minWorth < 0.0 || this.minWorth > 1.0)
            throw new ParseFailureException("Worthiness level must be between 0.0 and 1.0.");
    }

    @Override
    protected void validatePipeInput(TabbedLineReader inputStream) throws IOException {
        this.worthColIdx = inputStream.findField("worthiness");
        if (inputStream.findField("hammer") != 0)
            throw new IOException("Input does not look like a hammer file.  First column title must be \"hammer\".");
    }

    @Override
    protected void runPipeline(TabbedLineReader inputStream, PrintWriter writer) throws Exception {
        // Write the output header.
        writer.println(inputStream.header());
        // Loop through the input, copying to the output.
        int inCount = 0;
        int outCount = 0;
        for (var line : inputStream) {
            inCount++;
            if (line.getDouble(this.worthColIdx) >= this.minWorth) {
                outCount++;
                writer.println(line.toString());
            }
            if (log.isInfoEnabled() && inCount % 5000 == 0)
                log.info("{} lines read, {} written.", inCount, outCount);
        }
        log.info("{} lines read, {} written using minWorth = {}.", inCount, outCount, this.minWorth);
    }

}
