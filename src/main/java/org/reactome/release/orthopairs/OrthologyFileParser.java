package org.reactome.release.orthopairs;

import org.json.simple.JSONObject;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class OrthologyFileParser {
    private static Map<String,Map<String, Set<String>>> sourceTargetProteinHomologs = new HashMap<>();
    private static Map<String,Map<String, Set<String>>> targetGeneProteinMap = new HashMap<>();

    public static void parsePantherOrthologFiles(List<String> pantherFiles, String sourceMappingSpecies, JSONObject speciesJSONFile) throws IOException {

        System.out.println("Parsing homolog information from PANTHER files");
        // Panther uses different naming conventions for species, which needs to be mapped to Reactome's 4-letter species keys
        Set<String> pantherSpeciesNames = new HashSet<>();
        String sourceSpeciesPantherName = "";
        for (Object speciesKey : speciesJSONFile.keySet()) {
            JSONObject speciesJSON = (JSONObject) speciesJSONFile.get(speciesKey);
            if (!speciesKey.equals(sourceMappingSpecies)) {
                pantherSpeciesNames.add(speciesJSON.get("panther_name").toString());
            } else {
                sourceSpeciesPantherName = speciesJSON.get("panther_name").toString();
            }
        }

        if (!sourceSpeciesPantherName.equals("")) {
            // There are 2 Panther files at time of writing since Sus Scrofa (PIG) info is found in a seperate file.
            // The program iterates through both since for overlapping species it produces the same data structure. This is where the redundancy reduction of Sets comes in handy.
            for (String pantherFileTar : pantherFiles) {
                String extractedPantherFile = pantherFileTar.replace(".tar.gz", "");
                String line;
                BufferedReader br = new BufferedReader(new FileReader(extractedPantherFile));
                while ((line = br.readLine()) != null) {
                    // Sample line: HUMAN|HGNC=10663|UniProtKB=O60524	MOUSE|MGI=MGI=1918305|UniProtKB=Q8CCP0	LDO	Euarchontoglires	PTHR15239
                    // Tab-seperated sections: Human gene and protein info, Species gene and protein info, Ortholog type, last common ancestor for homolog (verify), pantherID
                    // For Ortholog type, we only want to look at the Least Divered Ortholog(LDO) and Ortholog(O) lines.
                    if (line.startsWith(sourceSpeciesPantherName)) {
                        String[] tabSplit = line.split("\t");

                        String sourceInfo = tabSplit[0];
                        String[] sourceSplit = sourceInfo.split("\\|");
                        String sourceGene = sourceSplit[1];
                        String sourceProtein = sourceSplit[2];

                        String targetInfo = tabSplit[1];
                        String[] targetSplit = targetInfo.split("\\|");
                        String targetSpecies = targetSplit[0];
                        String targetGene = targetSplit[1];
                        String targetProtein = targetSplit[2];

                        String orthologType = tabSplit[2];

                        // We don't look at lines that contain species that aren't in Reactome, or lines where the gene value starts with Gene (Gene|GeneID|Gene_Name|Gene_ORFName|Gene_OrderedLocusName)
                        // since these are often just names, not IDs. Additionally, we only want lines where orthologType is either an 'LDO' or 'O', and not a 'P', 'X', or 'LDX'.
                        if (pantherSpeciesNames.contains(targetSpecies) && !sourceGene.startsWith("Gene") && !targetGene.startsWith("Gene") && orthologType.contains("O")) {
                            sourceTargetProteinHomologs = MapId(targetSpecies, sourceProtein, targetProtein, sourceTargetProteinHomologs, orthologType);
                            targetGeneProteinMap = MapId(targetSpecies, targetGene, targetProtein, targetGeneProteinMap, orthologType);
                        }

                    }
                }
                br.close();
            }
        } else {
            System.out.println("Could not find source species in Species.json");
        }
    }

    private static Map<String, Map<String, Set<String>>> MapId(String targetSpecies, String keyEntity, String valueEntity, Map<String,Map<String,Set<String>>> entityMap, String orthologType) {

        // The 'null' conditionals are for the first time that level of the map is added. For example, if 'MOUSE' doesn't have an existing structure in the map,
        // the first condition will handle this.
        if (entityMap.get(targetSpecies) == null) {
            Set<String>  firstTargetEntityIdAdded = new HashSet<>(Arrays.asList(valueEntity));
            Map<String,Set<String>> firstSourceEntityMap = new HashMap<>();
            firstSourceEntityMap.put(keyEntity, firstTargetEntityIdAdded);
            entityMap.put(targetSpecies, firstSourceEntityMap);
        } else {
            if (entityMap.get(targetSpecies).get(keyEntity) == null) {
                Set<String> firstTargetEntityIdAdded = new HashSet<>(Arrays.asList(valueEntity));
                entityMap.get(targetSpecies).put(keyEntity, firstTargetEntityIdAdded);
            } else {
                entityMap.get(targetSpecies).get(keyEntity).add(valueEntity);
            }
        }

        // Lines with an orthologType equal to 'LDO' mean that we only want that value in the Set since its the Least Diverged Ortholog, meaning we have a
        // high degree of confidence in its homology. We remove all other values from the Set unless an 'LDO' value already exists. In this rare case,
        // we will keep multiple LDOs.
        Set<String> targetEntitys = entityMap.get(targetSpecies).get(keyEntity);
        if (orthologType.equals("LDO")) {
            if (!targetEntitys.contains("LDO")) {
                targetEntitys.clear();
                targetEntitys.add(valueEntity);
                targetEntitys.add("LDO");
                entityMap.get(targetSpecies).put(keyEntity, targetEntitys);
            }
        } else if (orthologType.equals("O") && targetEntitys.contains("LDO")) {
            targetEntitys.remove(valueEntity);
            entityMap.get(targetSpecies).put(keyEntity, targetEntitys);
        }
        return entityMap;
    }

    public static Map<String, Map<String, Set<String>>> getSourceAndTargetProteinHomologs() {
        return sourceTargetProteinHomologs;
    }

    public static Map<String, Map<String, Set<String>>> getTargetGeneProteinMap() {
        return targetGeneProteinMap;
    }
}