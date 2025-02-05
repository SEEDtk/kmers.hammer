/**
 *
 */
package org.theseed.sequence;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.util.Map;

import org.junit.jupiter.api.Test;
import org.theseed.sequence.SequenceDirectory.Member;

/**
 * @author Bruce Parrello
 *
 */
public class TestSequenceDir {

    /** flags used to verify that all members are found */
    final boolean membersFound[] = new boolean[8];
    /** map of member IDs to array indices */
    final static Map<String, Integer> indexMap = Map.of("1313.9775", 0, "171101.6", 1, "28181.21", 2,
            "1729695.3", 3, "1120961.3", 4, "1002870.3", 5, "1262806.3", 6, "2485170.3", 7);
    /** array of member names */
    final static String[] nameArray = new String[] { "Streptococcus pneumoniae strain 36_CMN",
            "Streptococcus pneumoniae R6", "Magnetovibrio blakemorei strain MV-1",
            "hammer.test.fa", "Ahrensia_kielensis_DSM_5890.fna", "Gaiella occulta strain F2-233",
            "Clostridium sp. CAG:433", "Granulicella sp. GAS466"};

    @Test
    public void testStream() {
        File inDir = new File("data", "seq_test");
        SequenceDirectory members = new SequenceDirectory(inDir);
        assertThat(members.size(), equalTo(8));
        members.parallelStream().parallel().forEach(x -> this.testMember(x));
        for (int i = 0; i < 8; i++)
            assertThat(Integer.toString(i), this.membersFound[i]);
    }

    /**
     * Verify that a member is correct.
     *
     * @param x		member to test
     */
    private void testMember(Member x) {
        // Get the member number.
        String memberId = x.getId();
        Integer idx = indexMap.get(x.getId());
        assertThat(memberId, idx, not(nullValue()));
        assertThat(memberId + " already visited.", ! this.membersFound[idx]);
        this.membersFound[idx] = true;
        assertThat(memberId, x.getName(), equalTo(nameArray[idx]));
    }

}
