/**
 *
 */
package org.theseed.proteins.hammer;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import org.junit.jupiter.api.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * @author Bruce Parrello
 *
 */
class HammerMapTest {

    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(HammerMapTest.class);

    @Test
    void test() {
        // Insure we have the normal word size before we check the masking math.
        int basis = Integer.SIZE;
        if (basis == 32) {
            assertThat(HammerMap.LOWER_BIT_MASK, equalTo(0x3FFFFFFF));
            assertThat(HammerMap.LOWER_BITS, equalTo(30));
            assertThat(HammerMap.LOWER_CHARS, equalTo(15));
        }
    }

}
