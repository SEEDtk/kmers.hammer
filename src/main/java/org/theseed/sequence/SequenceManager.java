/**
 *
 */
package org.theseed.sequence;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This object manages a list of sequences read from a FASTA file.  The actual iteration is managed by the
 * subclass, which either keeps all the sequences in memory, reads each one from a file every time, or
 * uses some other mechanism. Note that streaming is not supported because of the need to stay compatible
 * with sequential files.
 *
 * @author Bruce Parrello
 *
 */
public abstract class SequenceManager {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(SequenceManager.class);

    /**
     * This is the base class for the iterators we create.  It insures they are closed (if needed)
     */
    public abstract static class ResourceIter implements Iterator<Sequence>, AutoCloseable {

        @Override
        public abstract boolean hasNext();

        @Override
        public abstract Sequence next();

        @Override
        public abstract void close();

    }

    /**
     * Enumeration for the types of sequence managers supported.
     */
    public enum Type {
        /** cache the entire sequence set in memory */
        MEMORY {
            @Override
            public SequenceManager create(File file) throws IOException {
                return new CacheSequenceManager(file);
            }
        },
        /** read the sequences in from the file one at a time */
        FILE {
            @Override
            public SequenceManager create(File file) throws IOException {
                return new FileSequenceManager(file);
            }
        };

        /**
         * Create a sequence manager of this type using the specified file.
         *
         * @param file	FASTA file containing the sequences
         *
         * @return a sequence manager of the appropriate type
         *
         * @throws IOException
         */
        public abstract SequenceManager create(File file) throws IOException;
    }

    /**
     * Construct a sequence manager from a FASTA file.
     *
     * @param file		FASTA file to read for the sequences
     */
    public SequenceManager(File file) {
        // We don't do anything here, it simply insures that the subclass has to implement a special constructor
    }

    /**
     * @return a new, unique iteration object for this sequence set
     */
    public abstract ResourceIter newIterator();

}
