/**
 *
 */
package org.theseed.proteins.hammer;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;

import org.junit.jupiter.api.Test;

/**
 * @author Bruce Parrello
 *
 */
class TestSourMap {

    @Test
    void testSaveLoad() throws IOException {
        HammerSourMap testMap = HammerSourMap.create(new File("data", "tiny.hammers.tbl"));
        this.checkMap(testMap);
        File saveFile = new File("data", "tiny.sourmap.ser");
        testMap.save(saveFile);
        HammerSourMap loadMap = HammerSourMap.load(saveFile);
        this.checkMap(loadMap);
    }

    /**
     * Verify a SOUR map has all the correct role mappings.
     *
     * @param testMap
     */
    private void checkMap(HammerSourMap testMap) {
        assertThat("fig|2507161.3.peg.2267", testMap.getRole("fig|2507161.3.peg.2267"), equalTo("AlanTrnaSynt"));
        assertThat("fig|2692818.3.peg.3622", testMap.getRole("fig|2692818.3.peg.3622"), equalTo("AlanTrnaSynt"));
        assertThat("fig|2852100.3.peg.40", testMap.getRole("fig|2852100.3.peg.40"), equalTo("AlanTrnaSynt"));
        assertThat("fig|45361.3.peg.294", testMap.getRole("fig|45361.3.peg.294"), equalTo("AlanTrnaSynt"));
        assertThat("fig|412034.4.peg.370", testMap.getRole("fig|412034.4.peg.370"), equalTo("GlutTrnaSynt"));
        assertThat("fig|2173854.3.peg.450", testMap.getRole("fig|2173854.3.peg.450"), equalTo("IsolTrnaSynt"));
        assertThat("fig|2883571.10.peg.231", testMap.getRole("fig|2883571.10.peg.231"), equalTo("IsolTrnaSynt"));
        assertThat("fig|535091.3.peg.406", testMap.getRole("fig|535091.3.peg.406"), equalTo("LeucTrnaSynt"));
        assertThat("fig|997765.3.peg.149", testMap.getRole("fig|997765.3.peg.149"), equalTo("LeucTrnaSynt"));
        assertThat("fig|255719.3.peg.51", testMap.getRole("fig|255719.3.peg.51"), equalTo("NLThreSynt"));
        assertThat("fig|2852100.3.peg.373", testMap.getRole("fig|2852100.3.peg.373"), equalTo("SignRecoPartProt"));
        assertThat("fig|198804.5.peg.25", testMap.getRole("fig|198804.5.peg.25"), equalTo("SignRecoPartRece"));
        assertThat("fig|412034.4.peg.709", testMap.getRole("fig|412034.4.peg.709"), equalTo("SignRecoPartRece"));
        assertThat("fig|1258543.3.peg.122", testMap.getRole("fig|1258543.3.peg.122"), equalTo("ThreTrnaSynt"));
        assertThat("fig|192222.6.peg.134", testMap.getRole("fig|192222.6.peg.134"), equalTo("TranInitFact2n1"));
        assertThat("fig|2675773.3.peg.134", testMap.getRole("fig|2675773.3.peg.134"), equalTo("ValyTrnaSynt"));
        assertThat("fig|1925548.3.peg.589", testMap.getRole("fig|1925548.3.peg.589"), equalTo("ValyTrnaSynt"));
    }

}
