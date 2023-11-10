/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.IOException;

import org.kohsuke.args4j.Argument;
import org.theseed.basic.ParseFailureException;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.hammer.HammerSourMap;
import org.theseed.utils.BaseInputProcessor;

/**
 * This is a basic command that takes a hammer load file as input and outputs a SOUR map save file.
 * The resulting save file can be used in hammer analysis methods.
 *
 * The positional parameter is the name of the output save file.  The hammer load file should be on the
 * standard input.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -i	input hammer load file (if not STDIN)
 *
 * @author Bruce Parrello
 *
 */
public class HammerSourMapProcessor extends BaseInputProcessor {

    // COMMAND-LINE OPTIONS

    /** name of the output save file */
    @Argument(index = 0, metaVar = "saveFile", usage = "name of output file for the hammer SOUR map")
    private File saveFile;

    @Override
    protected void setReaderDefaults() {
    }

    @Override
    protected void validateReaderParms() throws IOException, ParseFailureException {
    }

    @Override
    protected void validateReaderInput(TabbedLineReader reader) throws IOException {
    }

    @Override
    protected void runReader(TabbedLineReader reader) throws Exception {
        HammerSourMap sourMap = HammerSourMap.create(reader);
        sourMap.save(saveFile);
    }

}
