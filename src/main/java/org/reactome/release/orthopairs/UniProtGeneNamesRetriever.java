package org.reactome.release.orthopairs;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import uk.ac.ebi.kraken.interfaces.uniprot.Gene;
import uk.ac.ebi.uniprot.dataservice.client.Client;
import uk.ac.ebi.uniprot.dataservice.client.QueryResult;
import uk.ac.ebi.uniprot.dataservice.client.QueryResultPage;
import uk.ac.ebi.uniprot.dataservice.client.ServiceFactory;
import uk.ac.ebi.uniprot.dataservice.client.exception.ServiceException;
import uk.ac.ebi.uniprot.dataservice.client.uniprot.UniProtComponent;
import uk.ac.ebi.uniprot.dataservice.client.uniprot.UniProtQueryBuilder;
import uk.ac.ebi.uniprot.dataservice.client.uniprot.UniProtService;
import uk.ac.ebi.uniprot.dataservice.query.Query;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.*;

public class UniProtGeneNamesRetriever {

    private static final Logger logger = LogManager.getLogger();
    private static final int MAX_NUMBER_OF_RETRIES = 5;
    private static final int MAX_UNIPROT_BATCH_QUERY_SIZE = 100;

    /**
     * Queries the UniProt mapping service through their Java API library. All Uniprot accession IDs are taken from the Panther
     * homolog files and used to query UniProt for the associated Gene Name. At time of writing (February 2020) this takes about
     * 3 hours to complete. The end result is a {targetSpecies}_gene_name_mapping.tsv file for each species.
     * @param speciesKey String - Shortened version of species name (eg: Bos taurus --> btau).
     * @param releaseNumber String - Used to create a directory where the produced files are stored.
     * @param pantherHomologMappings Map<String, Set<String>> - Species-specific homolog mappings.
     * @throws ServiceException - Thrown by UniProtService class if service is unavailable.
     * @throws InterruptedException - Thrown if the Sleep process that occurs after every batch query is interrupted.
     * @throws IOException - Thrown if unable to create/write to the mapping files.
     */
    public static void retrieveAndStoreGeneNameMappings(String speciesKey, String releaseNumber, Map<String, Set<String>> pantherHomologMappings) throws ServiceException, InterruptedException, IOException {

        // Get all gene names associated with all UniProt identifiers.
        Set<String> uniprotAccessionsToGeneNames = retrieveGeneNameMappings(pantherHomologMappings);
        // Store results to {targetSpecies}_gene_name_mapping.tsv file.
        storeGeneNameMappings(speciesKey, releaseNumber, uniprotAccessionsToGeneNames);
    }

    /**
     * Organizes all UniProt identifiers in the pantherMappings object into chunks of 250 identifiers. This is due to a
     * restriction in the size each query to UniProt can have.
     * @param pantherMappings Map<String, Set<String>> - Species-specific homolog mappings.
     * @return Set<String> - UniProt-Gene name mappings that will be written to file.
     * @throws ServiceException - Thrown by UniProtService class if service is unavailable.
     * @throws InterruptedException - Thrown if the Sleep process that occurs after every batch query is interrupted.
     */
    private static Set<String> retrieveGeneNameMappings(Map<String, Set<String>> pantherMappings) throws ServiceException, InterruptedException {
        // Get all UniProt identifiers and organize them into a List of Sets that contain 250 UniProt identifiers each.
        List<Set<String>> partitionedUniProtIds = getPartitionedUniProtIdentifiers(pantherMappings);
        //System.out.println(partitionedUniProtIds);
        // Query UniProt for gene names associated with UniProt identifiers.
        return retrieveGeneNamesFromUniProt(partitionedUniProtIds);

    }

    /**
     * Get all UniProt identifiers in pantherMappings and then partition them into chunks of 250 identifiers.
     * @param pantherMappings Map<String, Set<String>> - Species-specific homolog mappings.
     * @return <List><Set<String>> Partitioned UniProt identifiers.
     */
    private static List<Set<String>> getPartitionedUniProtIdentifiers(Map<String, Set<String>> pantherMappings) {

        // Get all UniProt identifiers in the Map, removing duplicates.
        Set<String> uniprotIds = getUniProtIdentifiers(pantherMappings);
        logger.info("Found " + uniprotIds.size() + " UniProt accessions");
        // Partition the UniProt identifiers into chunks of 250.
        return partitionUniProtIdentifiers(uniprotIds);
    }

    /**
     * Get all UniProt identifiers in pantherMappings.
     * @param pantherMappings Map<String, Set<String>> - Species-specific homolog mappings.
     * @return Set<String> - All UniProt identifiers in the mapping object.
     */
    public static Set<String> getUniProtIdentifiers(Map<String, Set<String>> pantherMappings) {
        Set<String> uniprotIds = new HashSet<>();
        for (String uniprotKey : pantherMappings.keySet()) {
            for (String uniprotValue : pantherMappings.get(uniprotKey)) {
                // Some identifiers aren't UniProt, and therefore can't be queried for Gene Names.
                if (uniprotValue.contains("UniProtKB")) {
                    // The values at this point are of the form "UniProtKB=123456", so we isolate the identifier.
                    uniprotIds.add(uniprotValue.split("=")[1]);
                }
            }
        }
        return uniprotIds;
    }

    /**
     * Partition all UniProt identifiers into 250-identifier chunks, and add each chunk to a List object.
     * @param uniprotIds Set<String> - All unique UniProt identifiers that were in the Panther mapping.
     * @return List<Set<String>> Partitioned UniProt identifiers.
     */
    public static List<Set<String>> partitionUniProtIdentifiers(Set<String> uniprotIds) {
        List<Set<String>> partitionedUniProtIds = new ArrayList<>();
        Set<String> partition = new HashSet<>();
        for (String uniprotId : uniprotIds) {
            // System.out.println(uniprotId);
            partition.add(uniprotId);
            // Once the partition has 250 identifiers, add it to the List and reset the partition variable.
            if (partition.size() == MAX_UNIPROT_BATCH_QUERY_SIZE) {
                partitionedUniProtIds.add(partition);
                partition = new HashSet<>();
            }
        }
        // The final identifiers at the end are added to the List here.
        partitionedUniProtIds.add(partition);
        return partitionedUniProtIds;
    }

    /**
     * Query UniProt via the Java API for gene names using UniProt accession IDs.
     * @param partitionedUniProtIds List<Set<String>> Partitioned UniProt identifiers.
     * @return Set<String> - Accession-to-Name mappings that will be stored to a file.
     * @throws ServiceException - Thrown by UniProtService class if service is unavailable.
     * @throws InterruptedException - Thrown if the Sleep process that occurs after every batch query is interrupted.
     */
    public static Set<String> retrieveGeneNamesFromUniProt(List<Set<String>> partitionedUniProtIds) throws ServiceException, InterruptedException {
        // Create the UniProtService.
        ServiceFactory serviceFactoryInstance = Client.getServiceFactoryInstance();
        UniProtService uniprotService = serviceFactoryInstance.getUniProtQueryService();
        uniprotService.start();
        Set<String> uniprotAccessionsToGeneNames = new HashSet<>();
        for (Set<String> uniprotIdentifierPartition : partitionedUniProtIds) {
            // Build UniProt API query from Set of 250 UniProt identifiers.
            //System.out.println(uniprotAccessionsToGeneNames);
            Query query = UniProtQueryBuilder.accessions(uniprotIdentifierPartition);
            // Perform UniProt API query to retrieve gene names associated with identifiers.

            QueryResult<UniProtComponent<Gene>> uniprotEntries = null;
            int currentQueryNum = 1;
            while (currentQueryNum < MAX_NUMBER_OF_RETRIES && uniprotEntries == null) {
                try {
                    uniprotEntries = uniprotService.getGenes(query);
                } catch (ServiceException e) {
                    logger.error(e);
                    currentQueryNum += 1;
                    Thread.sleep(5000);
                }
            }

            if (uniprotEntries == null) {
                throw new RuntimeException("Unable to get Uniprot Entries after " + MAX_NUMBER_OF_RETRIES + " tries");
            }


            //System.out.println("Number of hits: " + uniprotEntries.getNumberOfHits());
            QueryResultPage<UniProtComponent<Gene>> uniprotEntriesResultPage = uniprotEntries.getCurrentPage();
            int pageCount = 0;
            int totalHitCount = 0;
            while (totalHitCount < uniprotEntries.getNumberOfHits()) {
                //System.out.println(uniprotEntriesResultPage.pageSize());
                while (pageCount < uniprotEntriesResultPage.pageSize()) {
//                    if (count == 3409 || count == 3489 || count == 17775) {
//                        resultCount += 1;
//                        count += 1;
//                        continue;
//                    }
                    //System.out.println(totalHitCount);

                    UniProtComponent<Gene> geneObject;
                    try {
                        geneObject = uniprotEntriesResultPage.getResult(pageCount);
                        // System.out.println(geneObject.getAccession());
                    } catch (Exception e) {
                        pageCount += 1;
                        totalHitCount += 1;
                        continue;
                    }

                    //System.out.println(geneObject.getAccession() + ", Result: " + resultCount + ", Count: " + count);

                    pageCount += 1;
                    totalHitCount += 1;

                    //System.out.println(count);
                    //System.out.println(geneObject.getAccession());
                    if (!geneObject.getComponent().isEmpty()) {
                        // Iterate through all Gene components in the response.
                        for (Gene geneComponent : geneObject.getComponent()) {
                            // Tab-separate UniProt accession ID and its associated gene name, and then store these in the Set that will be returned.
                            uniprotAccessionsToGeneNames.add(geneObject.getAccession().toString() + "\t" + geneComponent.getGeneName().toString() + "\n");
                        }
                    }


                    if (totalHitCount % 1000 == 0) {
                        logger.info(totalHitCount + " UniProt identifiers have been queried for gene names");
                    }
                }
                uniprotEntriesResultPage = uniprotEntriesResultPage.fetchNextPage();
                //System.out.println("Page size: " + uniprotEntriesResultPage.pageSize());
                pageCount = 0;
            }

//            final AtomicInteger count1 = new AtomicInteger(0);
//uniprotEntries.forEachRemaining(geneObject -> {
//        count1.getAndIncrement();
//        System.out.println(geneObject.getAccession());
//        if (count1.get() % 1000 == 0) {
//            logger.info(count1.get() + " UniProt identifiers have been queried for gene names");
//        }
//    }
//);
//            List<UniProtComponent<Gene>> uniprotEntryList = IteratorUtils.toList(uniprotEntries);
//            for (UniProtComponent<Gene> uniprotEntry : uniprotEntryList) {
//                count++;
//                if (!uniprotEntry.getComponent().isEmpty()) {
//                    // Iterate through all Gene components in the response.
//                    for (Gene geneComponent : uniprotEntry.getComponent()) {
//                        // Tab-separate UniProt accession ID and its associated gene name, and then store these in the Set that will be returned.
//                        uniprotAccessionsToGeneNames.add(uniprotEntry.getAccession().toString() + "\t" + geneComponent.getGeneName().toString() + "\n");
//                    }
//
//
//                }
//
//                if (count % 1000 == 0) {
//                    logger.info(count + " UniProt identifiers have been queried for gene names");
//                }
//
//            }


            /*
            while (uniprotEntries.hasNext()) {
                count++;
                // Get Gene object returned from UniProt.
                try {
                    UniProtComponent<Gene> geneObject = uniprotEntries.next();
                    System.out.println(count);
                    System.out.println(geneObject.getAccession());
                    if (!geneObject.getComponent().isEmpty()) {
                        // Iterate through all Gene components in the response.
                        for (Gene geneComponent : geneObject.getComponent()) {
                            // Tab-separate UniProt accession ID and its associated gene name, and then store these in the Set that will be returned.
                            uniprotAccessionsToGeneNames.add(geneObject.getAccession().toString() + "\t" + geneComponent.getGeneName().toString() + "\n");
                        }
                    }
                } catch (Exception e) {
                    logger.error("Unable to get uniprot entry " + count);
                    uniprotEntries.next();
                    //System.exit(1);
                }

                if (count % 1000 == 0) {
                    logger.info(count + " UniProt identifiers have been queried for gene names");
                }
            }
            */
        }
        uniprotService.stop();

        return uniprotAccessionsToGeneNames;
    }

    /**
     * Write the UniProt-gene name mappings to the {targetSpecies}_gene_name_mapping.tsv file.
     * @param speciesKey String - Shortened version of species name.
     * @param releaseNumber String - Used to create a directory where the produced files are stored.
     * @param uniprotAccessionsToGeneNames Set<String> - Accession-to-Name mappings that will be written to a file.
     * @throws IOException - Thrown if unable to create/write to the mapping files.
     */
    private static void storeGeneNameMappings(String speciesKey, String releaseNumber, Set<String> uniprotAccessionsToGeneNames) throws IOException {

        Path uniprotAccessionsToGeneNamesFilePath = Paths.get(releaseNumber, speciesKey + "_gene_name_mapping.tsv");
        Files.deleteIfExists(uniprotAccessionsToGeneNamesFilePath);
        for (String uniprotAccessionsToGeneNamesLine : uniprotAccessionsToGeneNames) {
            Files.write(uniprotAccessionsToGeneNamesFilePath, uniprotAccessionsToGeneNamesLine.getBytes(), StandardOpenOption.CREATE, StandardOpenOption.APPEND);
        }
    }
}
