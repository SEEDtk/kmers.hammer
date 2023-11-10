/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.io.LineReader;
import org.theseed.io.TabbedLineReader;
import org.theseed.utils.BasePipeProcessor;

/**
 * This program takes as input a list of role IDs and a role definition file and outputs a smaller
 * role definition file containing only the specified role IDs.  Its primary use is to create the
 * SOUR definition file for "HammerProcessor".
 *
 * The positional parameter is the name of the roles.in.subsystems file (that is, the role definition
 * file).
 *
 * The list of role IDs will come in on the standard input.  It should be tab-delimited, with headers,
 * having the role IDs themselves in the first column.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -i	input file to use (if not STDIN)
 * -o	output file to use (if not STDOUT)
 *
 * @author Bruce Parrello
 *
 */
public class HammerDefineProcessor extends BasePipeProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(HammerDefineProcessor.class);
    /** set of roles to keep */
    private Set<String> goodRoles;
    /** expected number of roles */
    private static final int EXPECTED_ROLES = 100;

    // COMMAND-LINE OPTIONS

    /** role definition file */
    @Argument(index = 0, metaVar = "roles.in.subsystems", usage = "role definition file", required = true)
    private File roleFile;

    @Override
    protected void setPipeDefaults() {
    }

    @Override
    protected void validatePipeInput(TabbedLineReader inputStream) throws IOException {
        // Here we want to create our role set.
        this.goodRoles = new HashSet<String>(EXPECTED_ROLES);
        for (TabbedLineReader.Line line : inputStream)
            this.goodRoles.add(line.get(0));
        log.info("{} target role IDs found.", this.goodRoles.size());
    }

    @Override
    protected void validatePipeParms() throws IOException, ParseFailureException {
        // Verify the role definition file.
        if (! this.roleFile.canRead())
            throw new FileNotFoundException("Role definition file " + this.roleFile + " is not found or invalid.");
    }

    @Override
    protected void runPipeline(TabbedLineReader inputStream, PrintWriter writer) throws Exception {
        // Set up the counters and denote no roles have been found.
        int count = 0;
        int kept = 0;
        int skipped = 0;
        Set<String> foundRoles = new HashSet<String>(this.goodRoles.size() * 4 / 3);
        // Loop through the role definition file, keeping records we want.
        try (LineReader roleStream = new LineReader(this.roleFile)) {
            for (String line : roleStream) {
                count++;
                String role = StringUtils.substringBefore(line, "\t");
                if (! this.goodRoles.contains(role))
                    skipped++;
                else {
                    // Here we want to keep this line.
                    writer.println(line);
                    kept++;
                    foundRoles.add(role);
                }
            }
        }
        // Compute the roles not found.
        this.goodRoles.removeAll(foundRoles);
        if (this.goodRoles.size() > 0) {
            log.error("{} roles were not found.", this.goodRoles.size());
            log.info("Roles not found: {}", StringUtils.join(goodRoles, ", "));
        }
        log.info("{} lines read, {} kept, {} skipped.", count, kept, skipped);
    }

}
