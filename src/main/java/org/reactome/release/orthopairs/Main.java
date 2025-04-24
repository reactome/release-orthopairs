package org.reactome.release.orthopairs;

import java.io.*;
import java.nio.file.Paths;
import java.util.*;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

public class Main
{
    /* Orthopairs has been overhauled completely so that it now gets all Homology and Gene-Protein information from PANTHER (pantherdb.org) instead of Ensembl.
     *
     * In the past, Orthopairs was a 3-step process with the goal of getting Protein Homologs between the source species and all of Reactome's species:
     *  1) Download Gene-Protein information all species in Reactome and Protein-Gene information for the source species (Human in Reactome) from Ensembl's Biomart (An unstable step)
     *  2) Get Gene Homologs between source species and all species in Reactome from Ensembl Compara (A very long step)
     *  3) With the species Gene-Protein file from Biomart and the Gene Homolog file from Compara, map out a Protein homolog file
     *
     *  In this updated version, all of this information comes straight from 2 flat files provided by PANTHER. This makes things much quicker and stabler.
     *  The only twist to this though, is that not all of the species Gene IDs come from Ensembl, but rather there is some that come from organism-specific databases.
     *  The solution has been to use alternative ID mapping files from these organism-specific databases to convert any non-Ensembl IDs to their Ensembl equivalent.
     *
     *  At the end, we have the same two files for each species: The Gene-Protein mapping file specific to the species, and the Protein Homolog file.
     *  All Protein IDs are UniProt (instead of the mix of Ensembl and UniProt) and all Gene IDs match the format used by Ensembl.
     *
     *  These files are much smaller than the older version, since PANTHER annotates which homologs are the 'Least Diverged', meaning we can filter out many that add noise to the dataset.
     *
     */

    private static final Logger logger = LogManager.getLogger();

    public static void main( String[] args ) throws IOException, ParseException, InterruptedException {

        // If using an alternative source species, specify the 4-letter code as the second argument
        String pathToConfig = args.length > 0 ? args[0] : Paths.get("src", "main", "resources", "config.properties").toString();
        String sourceMappingSpecies = args.length > 1 ? args[1] : "hsap";
        Properties props = new Properties();
        props.load(new FileInputStream(pathToConfig));

        // Load config.properties
        String releaseNumber = props.getProperty("releaseNumber");
        String pathToSpeciesConfig = props.getProperty("pathToSpeciesConfig", "src/main/resources/Species.json");
        String pantherQfOFilename = props.getProperty("pantherQfOFilename", "QfO_Genome_Orthologs.tar.gz");
        String pantherHCOPFilename = props.getProperty("pantherHCOPFilename", "Orthologs_HCOP.tar.gz");

        if (releaseNumber.isEmpty()) {
            logger.fatal("Please populate config.properties file with releaseNumber");
            throw new IllegalStateException("No releaseNumber attribute in config.properties");
        }
        new File(releaseNumber).mkdir();

        logger.info("Starting Orthopairs file generation");
        List<String> pantherFiles = new ArrayList<String>(Arrays.asList(pantherQfOFilename, pantherHCOPFilename));

        JSONParser parser = new JSONParser();
        JSONObject speciesJSONFile = (JSONObject) parser.parse(new FileReader(pathToSpeciesConfig));

        // This method will produce two multi-level Maps, sourceTargetProteinHomologs and targetGeneProteinMap
        // sourceTargetProteinHomologs structure: {TargetSpecies-->{SourceProteinId-->[TargetHomologousProteinIds]}}
        // targetGeneProteinMap structure: {TargetSpecies-->{TargetGeneId-->[targetProteinIds]}}
        // The lower-level structure is a Set to reduce redundancy.
        OrthologyFileParser.parsePantherOrthologFiles(pantherFiles, sourceMappingSpecies, speciesJSONFile);
        Map<String,Map<String,Set<String>>> sourceTargetProteinHomologs = OrthologyFileParser.getSourceAndTargetProteinHomologs();
        Map<String,Map<String,Set<String>>> targetGeneProteinMap = OrthologyFileParser.getTargetGeneProteinMap();

        // Produces the protein homology, species gene-protein  and gene name mapping files
        for (Object speciesKey : speciesJSONFile.keySet()) {
            // No point in the source species mapping to itself
            if (!speciesKey.equals(sourceMappingSpecies)) {

                JSONObject speciesJSON = (JSONObject)speciesJSONFile.get(speciesKey);
                JSONArray speciesNames = (JSONArray) speciesJSON.get("name");
                logger.info("Attempting to create orthopairs files for " + speciesNames.get(0));
                String speciesPantherName = speciesJSON.get("panther_name").toString();
                // Produces the {sourceSpecies}_{targetspecies}_mapping.tsv file
                String sourceTargetProteinMappingFilename = releaseNumber + "/" + sourceMappingSpecies + "_" + speciesKey + "_mapping.tsv";
                Map<String,Set<String>> speciesProteinHomologs = sourceTargetProteinHomologs.get(speciesPantherName);
                OrthopairFileGenerator.createProteinHomologyFile(sourceTargetProteinMappingFilename, speciesProteinHomologs);
                // Produces the {targetSpecies}_gene_protein_mapping.tsv file
                String targetGeneProteinMappingFilename = releaseNumber + "/" + speciesKey + "_gene_protein_mapping.tsv";
                Map<String,Set<String>> speciesGeneProteinMap = targetGeneProteinMap.get(speciesPantherName);
                OrthopairFileGenerator.createSpeciesGeneProteinFile(speciesKey.toString(), targetGeneProteinMappingFilename, speciesJSON, speciesGeneProteinMap);
                // Queries UniProt API for gene names and creates the {targetSpecies}_gene_name_mapping.tsv file
                logger.info("Retrieving gene names for " + speciesNames.get(0) + " from UniProt");
                UniProtGeneNamesRetriever.retrieveAndStoreGeneNameMappings(speciesKey.toString(), releaseNumber, sourceTargetProteinHomologs.get(speciesPantherName));
            }
        }
        logger.info("Finished Orthopairs file generation");
    }
}
