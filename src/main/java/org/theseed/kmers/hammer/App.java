package org.theseed.kmers.hammer;

import java.util.Arrays;

import org.theseed.utils.BaseProcessor;

/**
 * This program processes DNA hammers.  The commands are as follows.
 *
 *  hammers		find all the hammers for a specified set of genomes
 */
public class App
{
    public static void main( String[] args )
    {
        // Get the control parameter.
        String command = args[0];
        String[] newArgs = Arrays.copyOfRange(args, 1, args.length);
        BaseProcessor processor;
        switch (command) {
        case "hammers" :
            processor = new HammerProcessor();
            break;
        default :
            throw new RuntimeException("Invalid command " + command + ".");
        }
        boolean ok = processor.parseCommand(newArgs);
        if (ok) {
            processor.run();
        }
    }
}
