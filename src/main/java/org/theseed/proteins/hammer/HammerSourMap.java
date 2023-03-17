/**
 *
 */
package org.theseed.proteins.hammer;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.HashMap;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.io.TabbedLineReader;

/**
 * This object maps a hammer FID to its associated role ID.  It can be created from a hammer load file
 * (such as produced by HammerFinderProcessor) or from a save file.
 *
 * @author Bruce Parrello
 *
 */
public class HammerSourMap implements Serializable {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(HammerSourMap.class);
    /** map from hammer feature IDs to role IDs */
    private HashMap<String, String> sourMap;
    /** object version ID for serialization */
    private static final long serialVersionUID = 2584000184555914262L;

    /**
     * Create a blank, empty sour map.
     */
    protected HammerSourMap() {
        this.sourMap = new HashMap<String, String>(1000);
    }

    /**
     * Build a sour map from a hammer load file.
     *
     * @param loadFile	name of the hammer load file
     *
     * @throws IOException
     */
    public static HammerSourMap create(File loadFile) throws IOException {
        log.info("Loading sour map from hammer file {}.", loadFile);
        HammerSourMap retVal;
        try (TabbedLineReader loadStream = new TabbedLineReader(loadFile)) {
            retVal = create(loadStream);
        }
        return retVal;
    }

    /**
     * Build a sour map from a hammer load file.  The hammer file is highly redundant, with each
     * feature ID appearing potentially hundreds of times.  We only update the hash the first
     * time the feature ID is encountered.
     *
     * @param loadStream	TabbedLineReader stream for the hammer load file
     *
     * @throws IOException
     */
    public static HammerSourMap create(TabbedLineReader loadStream) throws IOException {
        HammerSourMap retVal = new HammerSourMap();
        int fidIdx = loadStream.findField("fid");
        int roleIdx = loadStream.findField("role");
        long timer = System.currentTimeMillis();
        int inCount = 0;
        for (var line : loadStream) {
            String fid = line.get(fidIdx);
            retVal.sourMap.computeIfAbsent(fid, x -> line.get(roleIdx));
            inCount++;
            if (log.isInfoEnabled() && System.currentTimeMillis() - timer >= 5000) {
                timer = System.currentTimeMillis();
                log.info("{} hammers read.", inCount);
            }

        }
        log.info("{} hammer features found.", retVal.sourMap.size());
        return retVal;
    }

    /**
     * Load a saved sour map.
     *
     * @param saveFile	file to which the map was saved
     *
     * @throws IOException
     */
    public static HammerSourMap load(File saveFile) throws IOException {
        HammerSourMap retVal = null;
        try (FileInputStream inStream = new FileInputStream(saveFile)) {
            ObjectInputStream loadStream = new ObjectInputStream(inStream);
            retVal = (HammerSourMap) loadStream.readObject();
        } catch (ClassNotFoundException e) {
            throw new IOException(e.getMessage());
        }
        return retVal;
    }

    /**
     * Save this sour map to a file.
     *
     * @param saveFile	file to which the map should be saved
     *
     * @throws IOException
     */
    public void save(File saveFile) throws IOException {
        try (FileOutputStream outStream = new FileOutputStream(saveFile)) {
            ObjectOutputStream saveStream = new ObjectOutputStream(outStream);
            saveStream.writeObject(this);
            saveStream.close();
            log.info("Sour map saved to {}.", saveFile);
        }
    }

    /**
     * @return the role for a hammer FID, or an empty string if the FID is invalid
     *
     * @param hammerFid		feature ID of the hammer of interest
     */
    public String getRole(String hammerFid) {
        return this.sourMap.getOrDefault(hammerFid, "");
    }

    /**
     * @return the genome ID for a hammer FID
     *
     * @param hammerFid		feature ID of the hammer of interest
     */
    public String getRepId(String hammerFid) {
        return Feature.genomeOf(hammerFid);
    }

    /**
     * @return the number of hammer features in this map
     */
    public int size() {
        return this.sourMap.size();
    }

}
