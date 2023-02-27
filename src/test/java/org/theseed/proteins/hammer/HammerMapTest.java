/**
 *
 */
package org.theseed.proteins.hammer;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.junit.jupiter.api.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.io.TabbedLineReader;

/**
 * @author Bruce Parrello
 *
 */
class HammerMapTest {

    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(HammerMapTest.class);

    protected static class Thing {
        protected String fid;
        protected int count;

        public Thing(String id) {
            this.fid = id;
            this.count = 1;
        }
    }

    @Test
    void testBaseCounts() {
        Map<String, Integer> testMap = Map.of(
                "acgtacgtacgtac", 		 4,
                "acggtacggtaacggtacgt",	 7,
                "cccccccccc",			10,
                "acgtxACGTGAGC",		 4
        );
        for (var testEntry : testMap.entrySet()) {
            int count = HammerMap.commonBaseCount(testEntry.getKey());
            assertThat(testEntry.getKey(), count, equalTo(testEntry.getValue()));
        }

    }

    @Test
    void testMapParms() {
        // Insure we have the normal word size before we check the masking math.
        int basis = Integer.SIZE;
        if (basis == 32) {
            assertThat(HammerMap.LOWER_BIT_MASK, equalTo(0x3FFFFFFF));
            assertThat(HammerMap.LOWER_BITS, equalTo(30));
            assertThat(HammerMap.LOWER_CHARS, equalTo(15));
        }
    }

    @Test
    void testHammerMapSimple() {
        // Map the hammers to fid strings.
        var testMap = new HammerMap<String>(20);
        // Get a test hammer and its code.
        final String TEST_HAMMER = "acgtacgtacgtacgtacga";
        final long TEST_HAMMER_CODE = testMap.encode(TEST_HAMMER);
        // Build test hammers that go to the same sub-hash.
        final String TEST_HAMMER2 = "acgtacgtacgtacgtacgc";
        final String TEST_HAMMER3 = "acgtacgtacgtacgtacgg";
        final String TEST_HAMMER4 = "acgtacgtacgtacgtacgt";
        // Finally, a test hammer that chains to the first test-hammer
        final String TEST_HAMMERX = testMap.decode(TEST_HAMMER_CODE + 31);
        // Put all these hammers in the map.
        testMap.put(TEST_HAMMER, "test hammer");
        testMap.put(TEST_HAMMER2, "test hammer 2");
        testMap.put(TEST_HAMMER3, "test hammer 3");
        testMap.put(TEST_HAMMER4, "test hammer 4");
        testMap.put(TEST_HAMMERX, "test hammer X");
        // Check the size and the overload factor.
        assertThat(testMap.size(), equalTo(5L));
        assertThat(testMap.overloadFactor(), closeTo(1.25, 0.001));
        // Insure we can find the hammers.
        assertThat(testMap.get(TEST_HAMMER), equalTo("test hammer"));
        assertThat(testMap.get(TEST_HAMMER2), equalTo("test hammer 2"));
        assertThat(testMap.get(TEST_HAMMER3), equalTo("test hammer 3"));
        assertThat(testMap.get(TEST_HAMMER4), equalTo("test hammer 4"));
        assertThat(testMap.get(TEST_HAMMERX), equalTo("test hammer X"));
        assertThat(testMap.get("aaaaccccggggttttacgt"), nullValue());
        testMap.put("aaaaccccggggttttacgt", "alien hammer");
        assertThat(testMap.get("aaaaccccggggttttacgt"), equalTo("alien hammer"));
        int found = 0;
        for (var hammerInfo : testMap) {
            found++;
            String value = hammerInfo.getValue();
            switch (hammerInfo.getKey()) {
            case TEST_HAMMER :
                assertThat(value, equalTo("test hammer"));
                break;
            case TEST_HAMMER2 :
                assertThat(value, equalTo("test hammer 2"));
                break;
            case TEST_HAMMER3 :
                assertThat(value, equalTo("test hammer 3"));
                break;
            case TEST_HAMMER4 :
                assertThat(value, equalTo("test hammer 4"));
                break;
            case "aaaaccccggggttttacgt" :
                assertThat(value, equalTo("alien hammer"));
                break;
            default :
                assertThat(hammerInfo.getKey(), equalTo(TEST_HAMMERX));
                assertThat(value, equalTo("test hammer X"));
            }
        }
        assertThat(found, equalTo(6));
        final String ILLEGAL_HAMMER = "angtacgtacgtacgtacgt";
        assertThrows(IllegalArgumentException.class, () -> testMap.put(ILLEGAL_HAMMER, "illegal hammer"));
        assertThrows(IllegalArgumentException.class, () -> testMap.update(ILLEGAL_HAMMER, null, x -> x + " hammer"));
    }

    @Test
    void testHammerMapFile() throws IOException {
        var testMap = new HammerMap<Thing>(20);
        var checkMap = new HashMap<String, String>(30000);
        File inFile = new File("data", "hammers200.tbl");
        try (TabbedLineReader inStream = new TabbedLineReader(inFile)) {
            for (var line : inStream) {
                String hammer = line.get(0);
                String fid = line.get(1);
                checkMap.put(hammer, fid);
                testMap.put(hammer, new Thing(fid));
            }
        }
        log.info("Hammer map size is {}, overload factor is {}, load factor is {}.", testMap.size(), testMap.overloadFactor(),
                testMap.loadFactor());
        assertThat(testMap.size(), equalTo((long) checkMap.size()));
        for (var checkEntry : checkMap.entrySet()) {
            String hammer = checkEntry.getKey();
            Thing testThing = testMap.get(hammer);
            assertThat(testThing.fid, equalTo(checkEntry.getValue()));
        }
        boolean newKey = testMap.update("acgtacgtacgtaaccggtt", x -> x.count++, x -> new Thing(x));
        assertThat(newKey, equalTo(true));
        Thing thing = testMap.get("acgtacgtacgtaaccggtt");
        assertThat(thing.count, equalTo(1));
        Thing original = testMap.get("aaccttttaggtgtggaaaa");
        assertThat(original.count, equalTo(1));
        newKey = testMap.update("aaccttttaggtgtggaaaa", x -> x.count++, x -> new Thing(x));
        assertThat(newKey, equalTo(false));
        thing = testMap.get("aaccttttaggtgtggaaaa");
        assertThat(thing, sameInstance(original));
        assertThat(thing.count, equalTo(2));
    }


}
