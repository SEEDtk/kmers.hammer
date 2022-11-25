/**
 *
 */
package org.theseed.sequence;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

/**
 * @author Bruce Parrello
 *
 */
class TestSequenceManager {

    @Test
    void test() throws IOException {
        List<Sequence> list = new ArrayList<Sequence>(8);
        File contigs = new File("data", "contigs.fa");
        try (FastaInputStream inStream = new FastaInputStream(contigs)) {
            for (var seq : inStream)
                list.add(seq);
        }
        // We can use the list above to verify the results.
        SequenceManager mgr1 = SequenceManager.Type.FILE.create(contigs);
        SequenceManager mgr2 = SequenceManager.Type.MEMORY.create(contigs);
        // Get two iterators for each manager.  This insures we test parallelism.
        try (var iter1a = mgr1.newIterator(); var iter1b = mgr1.newIterator();
                var iter2a = mgr2.newIterator(); var iter2b = mgr2.newIterator()) {
            int count = 0;
            while (iter1a.hasNext()) {
                String message = String.format("After %d sequences.", count);
                assertThat(message, iter1b.hasNext(), equalTo(true));
                assertThat(message, iter2a.hasNext(), equalTo(true));
                assertThat(message, iter2b.hasNext(), equalTo(true));
                assertThat(message, iter1a.next(), equalTo(list.get(count)));
                assertThat(message, iter1b.next(), equalTo(list.get(count)));
                assertThat(message, iter2a.next(), equalTo(list.get(count)));
                assertThat(message, iter2b.next(), equalTo(list.get(count)));
                count++;
            }
            assertThat(iter1b.hasNext(), equalTo(false));
            assertThat(iter2a.hasNext(), equalTo(false));
            assertThat(iter2b.hasNext(), equalTo(false));
        }
    }

}
