/**
 *
 */
package org.theseed.proteins.hammer;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;
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
    /** nucleotide letters */
    private static final String[] NUCLEONS = new String[] { "a", "c", "g", "t" };

    protected static class Thing implements HammerMap.IScore {
        protected String fid;
        protected int count;
        protected boolean bad;

        public Thing(String id) {
            this.fid = id;
            this.count = 1;
            this.bad = false;
        }

        @Override
        public boolean isBadHammer() {
            return this.bad;
        }

        @Override
        public void setBadHammer() {
            this.bad = true;
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

    @Test
    public void testAnchorizing() throws IOException {
        var testMap = new HammerMap<Thing>(20);
        File inFile = new File("data", "hammers200.tbl");
        try (TabbedLineReader inStream = new TabbedLineReader(inFile)) {
            for (var line : inStream) {
                String hammer = line.get(0);
                String fid = line.get(1);
                testMap.put(hammer, new Thing(fid));
            }
        }
        int found = testMap.anchorize();
        log.info("{} non-anchors found.", found);
        var iter = testMap.iterator();
        int checkCount = 0;
        while (iter.hasNext()) {
            var hammerEntry = iter.next();
            checkCount++;
            String hammer = hammerEntry.getKey();
            if (hammerEntry.getValue().isBadHammer())
                this.checkAnchor(hammer, testMap);
        }
        log.info("{} hammers checked.", checkCount);
        // Now read the pure-anchor hammer map.
        testMap = new HammerMap<Thing>(20);
        inFile = new File("data", "anchors.tbl");
        try (TabbedLineReader inStream = new TabbedLineReader(inFile)) {
            for (var line : inStream) {
                String hammer = line.get(0);
                String fid = line.get(1);
                testMap.put(hammer, new Thing(fid));
            }
        }
        testMap.anchorize();
        iter = testMap.iterator();
        while (iter.hasNext()) {
            var hammerEntry = iter.next();
            String hammer = hammerEntry.getKey();
            assertThat(hammer, hammerEntry.getValue().isBadHammer(), equalTo(true));
        }
    }

    /**
     * Verify that the specified hammer has a close hammer in the map.
     *
     * @param hammer	hammer to test
     * @param testMap	map to test it against
     */
    private void checkAnchor(String hammer, HammerMap<Thing> testMap) {
        boolean found = false;
        final int n = hammer.length();
        for (int i = 0; ! found && i < n; i++) {
            String prefix = StringUtils.substring(hammer, 0, i);
            String suffix = StringUtils.substring(hammer, i+1);
            for (int j = 0; j < NUCLEONS.length; j++) {
                String newHammer = prefix + NUCLEONS[j] + suffix;
                if (! newHammer.equals(hammer)) {
                    Thing hammerData = testMap.get(newHammer);
                    if (hammerData != null)
                        found = true;
                }
            }
        }
        assertThat(hammer, found, equalTo(true));
    }


}
