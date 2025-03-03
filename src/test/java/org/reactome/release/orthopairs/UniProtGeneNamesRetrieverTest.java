package org.reactome.release.orthopairs;

import org.junit.Test;

import java.util.*;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

public class UniProtGeneNamesRetrieverTest {

    @Test
    public void testGetUniProtIdentifiers() {

        Set<String> testProteinIds = new HashSet<>(Arrays.asList("UniProtKB=12345"));
        Map<String, Set<String>> fakePantherMappings = new HashMap<>();
        fakePantherMappings.put("testKey", testProteinIds);

        Set<String> testUniProtIds = UniProtGeneNamesRetriever.getUniProtIdentifiers(fakePantherMappings);

        assertThat(testUniProtIds, hasSize(1));
        assertThat(testUniProtIds.contains("12345"), is(equalTo(true)));
    }

    @Test
    public void testGetUniProtIdentifiersReturnsOnlyUniProtValues() {
        Set<String> testProteinIds = new HashSet<>(Arrays.asList("UniProtKB=12345", "UniProtKB=67890", "NotARealKB=13579"));
        Map<String, Set<String>> fakePantherMappings = new HashMap<>();
        fakePantherMappings.put("testKey", testProteinIds);

        Set<String> testUniProtIds = UniProtGeneNamesRetriever.getUniProtIdentifiers(fakePantherMappings);

        assertThat(testUniProtIds, hasSize(2));
        assertThat(testUniProtIds.contains("13579"), is(equalTo(false)));
    }

    @Test
    public void testGetUniProtIdentifiersRemovesDuplicates() {
        Set<String> testProteinIds1 = new HashSet<>(Arrays.asList("UniProtKB=12345", "UniProtKB=67890"));
        Set<String> testProteinIds2 = new HashSet<>(Arrays.asList("UniProtKB=12345"));
        Map<String, Set<String>> fakePantherMappings = new HashMap<>();
        fakePantherMappings.put("testKey1", testProteinIds1);
        fakePantherMappings.put("testKey2", testProteinIds2);

        Set<String> testUniProtIds = UniProtGeneNamesRetriever.getUniProtIdentifiers(fakePantherMappings);

        assertThat(testUniProtIds, hasSize(2));
    }

    @Test
    public void testPartitionUniProtIdentifiers() {

        Set<String> testUniProtIds = new HashSet<>();
        for (int i = 0; i < 950; i++) {
            testUniProtIds.add("test" + i);
        }

        List<Set<String>> testPartitionedUniProtIds = UniProtGeneNamesRetriever.partitionUniProtIdentifiers(testUniProtIds);

        assertThat(testPartitionedUniProtIds, hasSize(4));
        // Since the total number of values added to the 'testUniProtIds' Set (950) does not cleanly divide
        // by 250 (max # of identifiers allowed per UniProt query), the last partition should be the remainder (200).
        assertThat(testPartitionedUniProtIds.get(3).size(), is(equalTo(200)));
    }

    @Test
    public void testPartitionUniProtIdentifiersRemovesDuplicates() {
        Set<String> testUniProtIds = new HashSet<>();
        for (int i = 0; i < 300; i++) {
            testUniProtIds.add("test");
        }

        List<Set<String>> testPartitionedUniProtIds = UniProtGeneNamesRetriever.partitionUniProtIdentifiers(testUniProtIds);

        assertThat(testPartitionedUniProtIds, hasSize(1));
        assertThat(testPartitionedUniProtIds.get(0).size(), is(equalTo(1)));
    }
}
