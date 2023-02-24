/**
 *
 */
package org.theseed.proteins.hammer;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;

import org.junit.jupiter.api.Test;
import org.theseed.io.TabbedLineReader;

/**
 * @author Bruce Parrello
 *
 */
class EncodingTest {

    @Test
    void testHammerCoding() throws IOException {
        File testFile = new File("data", "hammers200.tbl");
        try (TabbedLineReader testStream = new TabbedLineReader(testFile)) {
            long oldEncoded = -1;
            for (var line : testStream) {
                String hammerIn = line.get(0);
                long encoded = HammerMap.encode(hammerIn, 20);
                assertThat(hammerIn, encoded, not(equalTo(oldEncoded)));
                assertThat(hammerIn, encoded, greaterThanOrEqualTo(0l));
                String hammerOut = HammerMap.decode(encoded, 20);
                assertThat(hammerOut, equalTo(hammerIn));
                oldEncoded = encoded;
            }
        }
    }

}
