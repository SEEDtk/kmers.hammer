package org.theseed.kmers.hammer;

import java.util.Arrays;

import org.theseed.utils.BaseProcessor;

/**
 * This program processes DNA hammers.  The commands are as follows.
 *
 *  hammers		find all the hammers for a specified set of genomes
 *  closest		compute the closest genomes using the hammers
 *  define		create a hammer definition file from a roles.in.subsystems file
 *  dbLoad		load a hammer database from a flat file
 *  binTest		separate contigs into bins using hammers
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
        case "closest" :
            processor = new FindClosestProcessor();
            break;
        case "define" :
            processor = new HammerDefineProcessor();
            break;
        case "dbLoad" :
            processor = new HammerDbLoadProcessor();
            break;
        case "binTest" :
            processor = new BinTestProcessor();
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
