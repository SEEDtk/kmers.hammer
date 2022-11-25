/**
 *
 */
package org.theseed.sequence;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.UncheckedIOException;

/**
 * This sequence manager leaves the sequences on disk and reads through the file one record at a time.
 *
 * @author Bruce Parrello
 *
 */
public class FileSequenceManager extends SequenceManager {

    // FIELDS
    /** name of FASTA file */
    private File fastaFile;

    /**
     * Create an iterator through the file.
     *
     */
    protected class Iter extends SequenceManager.ResourceIter {

        /** fasta input stream */
        private FastaInputStream stream;

        /**
         * Create an iterator through the file.
         */
        public Iter() {
            try {
                this.stream = new FastaInputStream(FileSequenceManager.this.fastaFile);
            } catch (FileNotFoundException e) {
                throw new UncheckedIOException(e);
            }
        }

        @Override
        public boolean hasNext() {
            return this.stream.hasNext();
        }

        @Override
        public Sequence next() {
            return this.stream.next();
        }

        @Override
        public void close() {
            this.stream.close();
        }

    }

    /**
     * Construct a new file-based sequence manager.
     *
     * @param file	FASTA file containing the sequences
     *
     * @throws IOException
     */
    public FileSequenceManager(File file) throws IOException {
        super(file);
        // Verify the file exists.
        if (! file.canRead())
            throw new FileNotFoundException("FASTA file " + file + " not found or unreadable.");
        this.fastaFile = file;
    }

    @Override
    public ResourceIter newIterator() {
        return this.new Iter();
    }

}
