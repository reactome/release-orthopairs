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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.*;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class UniProtGeneNamesRetriever {

    private static final Logger logger = LogManager.getLogger();

    private static final String UNIPROT_ID_MAPPING_URL = "https://rest.uniprot.org/idmapping";
    private static final int MAX_NUMBER_OF_RETRIES = 5;
    private static final int UNIPROT_IDENTIFIER_BATCH_SIZE = 500;

    public enum RetrievalMethod {
        REST_API,
        JAVA_API
    }

    /**
     * Queries the UniProt mapping service through their REST API. All Uniprot accession IDs are taken from
     * the Panther homolog files and used to query UniProt for the associated Gene Name.  The end result is a
     * {targetSpecies}_gene_name_mapping.tsv file for each species.
     * @param speciesKey String - Shortened version of species name (eg: Bos taurus --> btau).
     * @param releaseNumber String - Used to create a directory where the produced files are stored.
     * @param pantherHomologMappings Map<String, Set<String>> - Species-specific homolog mappings.
     * @throws InterruptedException - Thrown if the Sleep process that occurs after every batch query is interrupted.
     * @throws IOException - Thrown if unable to create/write to the mapping files or if curl process to query
     * @throws ServiceException - Thrown by UniProtService class if service is unavailable.
     * REST API throws an exception
     */
    public static void retrieveAndStoreGeneNameMappings(
        String speciesKey, String releaseNumber, Map<String, Set<String>> pantherHomologMappings, RetrievalMethod retrievalMethod)
		throws IOException, InterruptedException, ServiceException {

        // Get all gene names associated with all UniProt identifiers.
        Set<String> uniprotAccessionsToGeneNames = retrieveGeneNameMappings(pantherHomologMappings, retrievalMethod);
        // Store results to {targetSpecies}_gene_name_mapping.tsv file.
        storeGeneNameMappings(speciesKey, releaseNumber, uniprotAccessionsToGeneNames);
    }

    /**
     * Organizes all UniProt identifiers in the pantherMappings object into chunks of identifiers. This is due to a
     * restriction in the size each query to UniProt can have.
     * @param pantherMappings Map<String, Set<String>> - Species-specific homolog mappings.
     * @return Set<String> - UniProt-Gene name mappings that will be written to file.
     * @throws IOException - Thrown if curl process to query REST API throws an exception
     * @throws InterruptedException - Thrown if the Sleep process that occurs after every batch query is interrupted.
     * @throws ServiceException - Thrown by UniProtService class if service is unavailable.
     */
    private static Set<String> retrieveGeneNameMappings(Map<String, Set<String>> pantherMappings, RetrievalMethod retrievalMethod)
		throws IOException, InterruptedException, ServiceException {
        // Get all UniProt identifiers and organize them into a List of Sets that contain 250 UniProt identifiers each.
        List<Set<String>> partitionedUniProtIds = getPartitionedUniProtIdentifiers(pantherMappings);
        //System.out.println(partitionedUniProtIds);
        // Query UniProt for gene names associated with UniProt identifiers.
        if (retrievalMethod.equals(RetrievalMethod.REST_API)) {
            return retrieveGeneNamesFromUniProtRESTAPI(partitionedUniProtIds);
        } else if (retrievalMethod.equals(RetrievalMethod.JAVA_API)) {
            return retrieveGeneNamesFromUniProtViaJavaAPI(partitionedUniProtIds);
        } else {
            throw new RuntimeException("UniProt Gene Name retrieval method unknown: " + retrievalMethod);
        }
    }

    /**
     * Get all UniProt identifiers in pantherMappings and then partition them into chunks of identifiers.
     * @param pantherMappings Map<String, Set<String>> - Species-specific homolog mappings.
     * @return <List><Set<String>> Partitioned UniProt identifiers.
     */
    private static List<Set<String>> getPartitionedUniProtIdentifiers(Map<String, Set<String>> pantherMappings) {
        // Get all UniProt identifiers in the Map, removing duplicates.
        Set<String> uniprotIds = getUniProtIdentifiers(pantherMappings);
        logger.info("Found " + uniprotIds.size() + " UniProt accessions");
        // Partition the UniProt identifiers into chunks.
        return partitionUniProtIdentifiers(uniprotIds);
    }

    /**
     * Get all UniProt identifiers in pantherMappings.
     * @param pantherMappings Map<String, Set<String>> - Species-specific homolog mappings.
     * @return Set<String> - All UniProt identifiers in the mapping object.
     */
    static Set<String> getUniProtIdentifiers(Map<String, Set<String>> pantherMappings) {
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
     * Partition all UniProt identifiers into identifier chunks, and add each chunk to a List object.
     * @param uniprotIds Set<String> - All unique UniProt identifiers that were in the Panther mapping.
     * @return List<Set<String>> Partitioned UniProt identifiers.
     */
    static List<Set<String>> partitionUniProtIdentifiers(Set<String> uniprotIds) {
        List<Set<String>> partitionedUniProtIds = new ArrayList<>();
        Set<String> partition = new HashSet<>();
        for (String uniprotId : uniprotIds) {
            partition.add(uniprotId);
            // Once the partition has reached the batch size, add it to the List and reset the partition variable.
            if (partition.size() == UNIPROT_IDENTIFIER_BATCH_SIZE) {
                partitionedUniProtIds.add(partition);
                partition = new HashSet<>();
            }
        }
        // The final identifiers at the end are added to the List here.
        partitionedUniProtIds.add(partition);
        return partitionedUniProtIds;
    }

    /**
     * Query UniProt via the REST API for gene names using UniProt accession IDs.
     * @param partitionedUniProtIds List<Set<String>> Partitioned UniProt identifiers.
     * @return Set<String> - Accession-to-Name mappings that will be stored to a file.
     * @throws IOException - Thrown if curl process to query REST API throws an exception
     * @throws InterruptedException - Thrown if the Sleep process that occurs after every batch query is interrupted.
     */
    private static Set<String> retrieveGeneNamesFromUniProtRESTAPI(List<Set<String>> partitionedUniProtIds) throws
		IOException, InterruptedException {

        Set<String> uniprotAccessionsToGeneNames = new HashSet<>();

        for (Set<String> uniProtIdentifiersBatch : partitionedUniProtIds) {

            String jobId = submitJob(uniProtIdentifiersBatch);
            while (!jobFinished(jobId)) {
                TimeUnit.SECONDS.sleep(1);
            }

            Map<String, String> results = fetchResults(jobId);
            for (Map.Entry<String,String> result : results.entrySet()) {
                String accession = result.getKey();
                String geneName = result.getValue();

                String uniprotAccessionToGeneName = accession + "\t" + geneName + "\n";
                uniprotAccessionsToGeneNames.add(uniprotAccessionToGeneName);
            }

            for (String uniProtIdentifier : uniProtIdentifiersBatch) {
                if (!results.containsKey(uniProtIdentifier)) {
                    Files.write(
                        Paths.get("not_mapped.tsv"),
                        uniProtIdentifier.concat(System.lineSeparator()).getBytes(),
                        StandardOpenOption.CREATE, StandardOpenOption.APPEND
                    );
                }
            }
        }

        return uniprotAccessionsToGeneNames;
    }

    /**
     * Query UniProt via the Java API for gene names using UniProt accession IDs.
     * @param partitionedUniProtIds List<Set<String>> Partitioned UniProt identifiers.
     * @return Set<String> - Accession-to-Name mappings that will be stored to a file.
     * @throws ServiceException - Thrown by UniProtService class if service is unavailable.
     * @throws InterruptedException - Thrown if the Sleep process that occurs after every batch query is interrupted.
     */
    private static Set<String> retrieveGeneNamesFromUniProtViaJavaAPI(List<Set<String>> partitionedUniProtIds) throws ServiceException, InterruptedException {
        // Create the UniProtService.
        ServiceFactory serviceFactoryInstance = Client.getServiceFactoryInstance();
        UniProtService uniprotService = serviceFactoryInstance.getUniProtQueryService();
        uniprotService.start();
        Set<String> uniprotAccessionsToGeneNames = new HashSet<>();
        for (Set<String> uniprotIdentifierPartition : partitionedUniProtIds) {
            // Build UniProt API query from Set of 250 UniProt identifiers.
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


            QueryResultPage<UniProtComponent<Gene>> uniprotEntriesResultPage = uniprotEntries.getCurrentPage();
            int pageCount = 0;
            int totalHitCount = 0;
            while (totalHitCount < uniprotEntries.getNumberOfHits()) {
                while (pageCount < uniprotEntriesResultPage.pageSize()) {

                    UniProtComponent<Gene> geneObject;
                    try {
                        geneObject = uniprotEntriesResultPage.getResult(pageCount);
                    } catch (Exception e) {
                        pageCount += 1;
                        totalHitCount += 1;
                        continue;
                    }

                    pageCount += 1;
                    totalHitCount += 1;

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
                pageCount = 0;
            }

        }
        uniprotService.stop();

        return uniprotAccessionsToGeneNames;
    }

    private static String submitJob(Set<String> accessions) throws IOException {
        StringBuilder curlQueryBuilder = new StringBuilder();
        curlQueryBuilder.append("curl --request POST ");
        curlQueryBuilder.append(UNIPROT_ID_MAPPING_URL + "/run ");
        curlQueryBuilder.append("--form ids=\"" + String.join( ",",accessions) + "\"" + " ");
        curlQueryBuilder.append("--form from=\"UniProtKB_AC-ID\" ");
        curlQueryBuilder.append("--form to=\"UniProtKB\" ");

        Process process = Runtime.getRuntime().exec(curlQueryBuilder.toString());
        BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(process.getInputStream()));

        return bufferedReader
            .lines()
            .filter(line -> line.contains("jobId"))
            .map(line -> {
                Matcher jobIdMatcher = Pattern.compile("\"jobId\":\"(.*)\"").matcher(line);
                if (jobIdMatcher.find()) {
                    return jobIdMatcher.group(1);
                } else {
                    throw new RuntimeException("Could not get job id from " + line);
                }
            })
            .findFirst()
            .orElseThrow(() -> new RuntimeException("Curl could not get job id"));
    }

    private static boolean jobFinished(String jobID) throws IOException {
        StringBuilder curlQueryBuilder = new StringBuilder();
        curlQueryBuilder.append("curl -s ");
        curlQueryBuilder.append(UNIPROT_ID_MAPPING_URL + "/status/" + jobID);

        Process process = Runtime.getRuntime().exec(curlQueryBuilder.toString());
        BufferedReader jobStatusWebSource = new BufferedReader(new InputStreamReader(process.getInputStream()));

        return jobStatusWebSource
            .lines()
            .anyMatch(
                line -> line.contains("{\"jobStatus\":\"FINISHED\"}")
            );
    }

    private static Map<String, String> fetchResults(String jobId) throws IOException {
        Map<String, String> accessionToGeneName = new HashMap<>();

        String jobResultsURL = UNIPROT_ID_MAPPING_URL + "/uniprotkb/results/stream/" + jobId +
            "?fields=gene_primary&format=tsv";

        Process process = Runtime.getRuntime().exec("curl -s " + jobResultsURL);
        BufferedReader jobResultsWebSource = new BufferedReader(new InputStreamReader(process.getInputStream()));

        jobResultsWebSource.readLine(); // Skip header

        String resultsLine;
        while ((resultsLine = jobResultsWebSource.readLine()) != null) {
            String[] columns = resultsLine.split("\t");
            if (columns.length == 2) {
                String accession = columns[0];
                String geneName = columns[1];

                accessionToGeneName.put(accession, geneName);
            }
        }
        return accessionToGeneName;
    }

    /**
     * Write the UniProt-gene name mappings to the {targetSpecies}_gene_name_mapping.tsv file.
     * @param speciesKey String - Shortened version of species name.
     * @param releaseNumber String - Used to create a directory where the produced files are stored.
     * @param uniprotAccessionsToGeneNames Set<String> - Accession-to-Name mappings that will be written to a file.
     * @throws IOException - Thrown if unable to create/write to the mapping files.
     */
    private static void storeGeneNameMappings(
        String speciesKey, String releaseNumber, Set<String> uniprotAccessionsToGeneNames) throws IOException {

        Path uniprotAccessionsToGeneNamesFilePath =
            Paths.get(releaseNumber, speciesKey + "_gene_name_mapping.tsv");

        Files.deleteIfExists(uniprotAccessionsToGeneNamesFilePath);

        for (String uniprotAccessionsToGeneNamesLine : uniprotAccessionsToGeneNames) {
            Files.write(
                uniprotAccessionsToGeneNamesFilePath,
                uniprotAccessionsToGeneNamesLine.getBytes(),
                StandardOpenOption.CREATE, StandardOpenOption.APPEND
            );
        }
    }
}
