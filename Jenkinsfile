// This Jenkinsfile is used by Jenkins to run the Orthopairs step of Reactome's release.
// This step downloads Orthology data from pantherdb.org and generates a set of orthology mapping files 
// for each species projected to. It also finds the gene names associated with each gene identifier via UniProt. 

import org.reactome.release.jenkins.utilities.Utilities

// Shared library maintained at 'release-jenkins-utils' repository.
def utils = new Utilities()

pipeline{
	agent any

	stages{
		// This stage checks that an upstream step, ConfirmReleaseConfigs, was run successfully.
		stage('Check ConfirmReleaseConfigs build succeeded'){
			steps{
				script{
					utils.checkUpstreamBuildsSucceeded("ConfirmReleaseConfigs")
				}
			}
		}
		// Download two orthology files, Orthologs_HCOP.tar.gz and QfO_Genome_Orthologs.tar.gz, generated by PANTHER used to create the Reactome orthology mappings.
		stage('Setup: Download ortholog files from PANTHER'){
			steps{
				script{
					def hcopFilename = "Orthologs_HCOP.tar.gz";
					def qfoFilename = "QfO_Genome_Orthologs.tar.gz";
					def pantherReleaseURL = "ftp://ftp.pantherdb.org/ortholog/current_release"
					sh "wget -q ${pantherReleaseURL}/${hcopFilename}"
					sh "tar -xvf ${hcopFilename}"
					sh "wget -q ${pantherReleaseURL}/${qfoFilename}"
					sh "tar -xvf ${qfoFilename}"
				}
			}
		}
		// Download files that contain mapping links between UniProt identifiers and identifiers specific to model organism dbs.
		stage('Setup: Download alternate ID mapping files from model organism databases'){
			steps{
				script{
					sh "wget -q -O mmus_alternate_ids.txt http://www.informatics.jax.org/downloads/reports/HGNC_AllianceHomology.rpt"
					sh "wget -q -O rnor_alternate_ids.txt ftp://ftp.rgd.mcw.edu/pub/data_release/GENES_RAT.txt"
					sh "wget -q -O xtro_alternate_ids.txt ftp://ftp.xenbase.org/pub/GenePageReports/GenePageEnsemblModelMapping.txt"
					sh "wget -q -O drer_alternate_ids.txt https://zfin.org/downloads/ensembl_1_to_1.txt"
				}
			}
		}	
		// This stage builds the jar file using maven.
		stage('Setup: Build jar file') {
			steps {
				script {
					utils.buildJarFile()
				}
			}
		}
		// This stage executes the Orthopairs jar file, producing all Orthopairs files used by Orthoinference.
		stage('Main: Generate Orthopairs files') {
			steps {
				script {
					// The credentials used here are a config file uploaded to Jenkins.
					withCredentials([file(credentialsId: 'Config', variable: 'ConfigFile')]) {
						sh "java -jar target/orthopairs-*-jar-with-dependencies.jar $ConfigFile"
					}
				}
			}
		}
		// This stage compares the line counts of the orthopairs files generated between the current and previous release.
		// An intelligible output should be visible at the console logs for the build.
		stage('Post: Orthopairs file line counts') {
		    steps {
		        script {
		            def releaseVersion = utils.getReleaseVersion()
		            def previousReleaseVersion = utils.getPreviousReleaseVersion()
		            def currentDir = pwd()
		            sh "mkdir -p ${previousReleaseVersion}/"
		            sh "aws s3 --recursive --no-progress cp s3://reactome/private/releases/${previousReleaseVersion}/orthopairs/data/orthopairs/ ${previousReleaseVersion}/"
		            sh "gunzip ${previousReleaseVersion}/*"
		            utils.outputLineCountsOfFilesBetweenFolders("$releaseVersion", "$previousReleaseVersion", "$currentDir")
		            sh "rm -r ${previousReleaseVersion}"
		        }
		    }
		}
		// Logs and data files generated by this step are archived in the Reactome S3 bucket.
		// All files are then deleted on server.
		stage('Post: Archive Outputs'){
			steps{
				script{
					def releaseVersion = utils.getReleaseVersion()
					sh "mkdir -p orthopairs/"
					sh "mv ${releaseVersion}/* orthopairs/"
					def dataFiles = ["orthopairs", "*alternate_ids.txt", "*.gz"]
					// Log files are automatically output to a 'logs' folder, so nothing needs to be specified here.
					def logFiles = []
					def foldersToDelete = []
					utils.cleanUpAndArchiveBuildFiles("orthopairs", dataFiles, logFiles, foldersToDelete)
					sh "rm *Orthologs*"
				}
			}
		}
	}
}
