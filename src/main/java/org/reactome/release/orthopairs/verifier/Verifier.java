package org.reactome.release.orthopairs.verifier;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;

import org.reactome.release.verifier.Results;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

import static org.reactome.release.verifier.CountUtils.greaterThanOrEqualTo5PercentDrop;
import static org.reactome.release.verifier.FileUtils.*;

/**
 * @author Joel Weiser (joel.weiser@oicr.on.ca)
 * Created 8/15/2024
 */
public class Verifier {
    @Parameter(names = "-release", description = "Current Reactome Release number", required = true)
    private int releaseNumber;

    public static void main(String[] args) throws IOException {
        Verifier verifier = new Verifier();
        JCommander.newBuilder()
            .addObject(verifier)
            .build()
            .parse(args);

        verifier.run();
    }

    private void run() throws IOException {
        Results results = compareCurrentAndPreviousOrthopairFiles();

        results.report();
        if (results.hasErrors()) {
            results.reportErrors();
            System.exit(1);
        }
    }

    private Results compareCurrentAndPreviousOrthopairFiles() throws IOException {
        Results overallResults = new Results();

        List<Path> previousOrthopairFilePaths = getPreviousOrthopairFilePaths();
        for (Path currentOrthopairFilePath : getCurrentOrthopairFilePaths()) {
            Path previousOrthopairFilePath =
                getAnalgousPreviousOrthopairFile(currentOrthopairFilePath, previousOrthopairFilePaths);

            if (previousOrthopairFilePath != null) {
                Results lineCountResults = checkLineCount(currentOrthopairFilePath, previousOrthopairFilePath);
                Results fileSizeResults = checkFileSize(currentOrthopairFilePath, previousOrthopairFilePath);

                overallResults.addInfoMessages(lineCountResults.getInfoMessages());
                overallResults.addInfoMessages(fileSizeResults.getInfoMessages());
                overallResults.addErrorMessages(lineCountResults.getErrorMessages());
                overallResults.addErrorMessages(fileSizeResults.getErrorMessages());
            }
        }

        return overallResults;
    }

    private List<Path> getCurrentOrthopairFilePaths() throws IOException {
        return Files.list(Paths.get(String.valueOf(releaseNumber))).collect(Collectors.toList());
    }

    private List<Path> getPreviousOrthopairFilePaths() throws IOException {
        String folderName = downloadPreviousOrthopairFiles();
        Files.list(Paths.get(folderName)).filter(file -> file.getFileName().toString().endsWith(".gz")).forEach(file -> {
            try {
                gunzipFile(file);
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        });
        return Files.list(Paths.get(folderName)).filter(file -> file.getFileName().toString().endsWith(".tsv")).collect(Collectors.toList());
    }

    private Path getAnalgousPreviousOrthopairFile(
        Path currentOrthopairFilePath, List<Path> previousOrthopairFilePaths) {

        return previousOrthopairFilePaths
            .stream()
            .filter(previousOrthopairFilePath ->
                previousOrthopairFilePath.getFileName().toString()
                    .equals(currentOrthopairFilePath.getFileName().toString())
            )
            .findFirst()
            .orElse(null);
    }

    private Results checkLineCount(Path currentOrthopairFilePath, Path previousOrthopairFilePath)
        throws IOException {

        Results lineCountResults = new Results();

        int currentLineCount = (int) Files.lines(currentOrthopairFilePath).count();
        int previousLineCount = (int) Files.lines(previousOrthopairFilePath).count();

        if (greaterThanOrEqualTo5PercentDrop(currentLineCount, previousLineCount)) {
            lineCountResults.addErrorMessage(
                String.format("%s has significantly fewer lines than %s (%d vs. %d) - Difference: %d lines",
                    currentOrthopairFilePath, previousOrthopairFilePath, currentLineCount, previousLineCount,
                    (previousLineCount - currentLineCount))
            );
        } else {
            lineCountResults.addInfoMessage(
                String.format("%s (%d lines) vs. %s (%d lines) - Difference: %d lines",
                    currentOrthopairFilePath, currentLineCount, previousOrthopairFilePath, previousLineCount,
                    (currentLineCount - previousLineCount))
            );
        }

        return lineCountResults;
    }

    private Results checkFileSize(Path currentOrthopairFilePath, Path previousOrthopairFilePath) throws IOException {

        Results fileSizeResults = new Results();

        int currentFileSize = (int) Files.size(currentOrthopairFilePath);
        int previousFileSize = (int) Files.size(previousOrthopairFilePath);

        if (greaterThanOrEqualTo5PercentDrop(currentFileSize, previousFileSize)) {
            fileSizeResults.addErrorMessage(
                String.format("%s has significantly smaller size than %s (%d bytes vs. %d bytes) - Difference: %d bytes",
                    currentOrthopairFilePath, previousOrthopairFilePath, currentFileSize, previousFileSize,
                    (previousFileSize - currentFileSize))
            );
        } else {
            fileSizeResults.addInfoMessage(
                String.format("%s (%d bytes) vs. %s (%d bytes) - Difference: %d bytes",
                    currentOrthopairFilePath, currentFileSize, previousOrthopairFilePath, previousFileSize,
                    (currentFileSize - previousFileSize))
            );
        }

        return fileSizeResults;
    }


    private String downloadPreviousOrthopairFiles() throws IOException {
        deleteDirectory(getPreviousReleaseNumberDirectory());
        downloadFolderFromS3("reactome", getS3FolderPath());
        return renameFolderToPreviousReleaseVersionNumber().toString();
    }

    private Path renameFolderToPreviousReleaseVersionNumber() throws IOException {
        Files.move(Paths.get(getS3FolderPath()).getFileName(), getPreviousReleaseNumberDirectory());
        return getPreviousReleaseNumberDirectory();
    }

    private String getS3FolderPath() {
        return String.format("private/releases/%d/orthopairs/data/orthopairs/", getPreviousReleaseNumber());
    }

    private Path getPreviousReleaseNumberDirectory() {
        return Paths.get(String.valueOf(getPreviousReleaseNumber()));
    }

    private int getPreviousReleaseNumber() {
        return this.releaseNumber - 1;
    }
}
