/**
 *
 */
package org.theseed.sequence.seeds.filters;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

import org.junit.jupiter.api.Test;
import org.theseed.genome.Feature;
import org.theseed.sequence.FastaInputStream;

/**
 * @author Bruce Parrello
 *
 */
class HammerFilterTests implements HammerDupFilter.IParms, HammerFeatureFilter.IParms {

    @Test
    void testDupFilters() throws IOException {

        File testFile = new File("data", "hammerFilters.fa");
        // Set up the testing filters.
        HammerDupFilter dupNone = HammerDupFilter.Type.NONE.create(this);
        HammerDupFilter dupKeep = HammerDupFilter.Type.KEEP.create(this);
        HammerDupFilter dupLongest = HammerDupFilter.Type.LONGEST.create(this);
        HammerFeatureFilter pegSimple = HammerFeatureFilter.Type.SIMPLE.create(this);
        HammerFeatureFilter pegAll = HammerFeatureFilter.Type.ALL.create(this);
        // We will set up the dups in here.
        Map<String, Map<String, String>> genomeFidMapMap = new HashMap<String, Map<String, String>>();
        try (var inStream = new FastaInputStream(testFile)) {
            for (var seq : inStream) {
                String fid = seq.getLabel();
                // Test the peg filters.
                boolean pegSimpleFlag = pegSimple.check(seq);
                boolean pegAllFlag = pegAll.check(seq);
                assertThat(fid, pegAllFlag, equalTo(true));
                switch (fid) {
                case "fig|2813572.3.peg.1735" :
                case "fig|1733.7982.peg.5278" :
                    assertThat(fid, pegSimpleFlag, equalTo(false));
                    break;
                default :
                    assertThat(fid, pegSimpleFlag, equalTo(true));
                    break;
                }
                // Put this sequence in the fid map.
                String gid = Feature.genomeOf(fid);
                Map<String, String> fidMap = genomeFidMapMap.computeIfAbsent(gid, x -> new TreeMap<String, String>());
                fidMap.put(fid, seq.getSequence());
            }
        }
        // Now process the fidMap.  There are two interesting cases-- 1733.7982 and 1733.8144.
        for (var genomeEntry : genomeFidMapMap.entrySet()) {
            String gid = genomeEntry.getKey();
            Map<String, String> fidMap = genomeEntry.getValue();
            var fidMapNone = dupNone.filter(fidMap);
            var fidMapKeep = dupKeep.filter(fidMap);
            var fidMapLongest = dupLongest.filter(fidMap);
            switch (gid) {
            case "1733.7982" :
                assertThat(fidMapNone.size(), equalTo(0));
                assertThat(fidMapKeep.size(), equalTo(2));
                for (String fid2 : fidMapKeep.keySet())
                    assertThat(fid2, fidMapKeep.get(fid2), equalTo(fidMap.get(fid2)));
                assertThat(fidMapLongest.size(), equalTo(1));
                assertThat(fidMapLongest.containsKey("fig|1733.7982.peg.4522"), equalTo(true));
                assertThat(fidMapLongest.get("fig|1733.7982.peg.4522"), equalTo(fidMap.get("fig|1733.7982.peg.4522")));
                break;
            case "1733.8144" :
                assertThat(fidMapNone.size(), equalTo(0));
                assertThat(fidMapKeep.size(), equalTo(2));
                for (String fid2 : fidMapKeep.keySet())
                    assertThat(fid2, fidMapKeep.get(fid2), equalTo(fidMap.get(fid2)));
                assertThat(fidMapLongest.size(), equalTo(1));
                assertThat(fidMapLongest.containsKey("fig|1733.8144.peg.1814"), equalTo(true));
                assertThat(fidMapLongest.get("fig|1733.8144.peg.1814"), equalTo(fidMap.get("fig|1733.8144.peg.1814")));
                break;
            default :
                assertThat(fidMapNone.size(), equalTo(1));
                assertThat(fidMapKeep.size(), equalTo(1));
                assertThat(fidMapLongest.size(), equalTo(1));
                for (String fid2 : fidMapNone.keySet()) {
                    assertThat(fid2, fidMapNone.get(fid2), equalTo(fidMap.get(fid2)));
                    assertThat(fid2, fidMapKeep.get(fid2), equalTo(fidMap.get(fid2)));
                    assertThat(fid2, fidMapLongest.get(fid2), equalTo(fidMap.get(fid2)));
                }
            }

        }

    }

}
