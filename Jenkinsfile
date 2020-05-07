import groovy.json.JsonSlurper
// This Jenkinsfile is used by Jenkins to run the Orthopairs step of Reactome's release.
// This step downloads Orthology data from pantherdb.org and generates a set of orthology mapping files 
// for each species projected to. It also finds the gene names associated with each gene identifier via UniProt. 
// It requires that the 'ConfirmReleaseConfigs' step has been run successfully before it can be run.
pipeline{
	agent any

	stages{
		// This stage checks that an upstream project, ConfirmReleaseConfig, was run successfully for its last build.
		stage('Check ConfirmReleaseConfig build succeeded'){
			steps{
				script{
					def currentRelease = (pwd() =~ /(\d+)\//)[0][1];
					// This queries the Jenkins API to confirm that the most recent build of 'ConfirmReleaseConfigs' was successful.
					def configStatusUrl = httpRequest authentication: 'jenkinsKey', validResponseCodes: "${env.VALID_RESPONSE_CODES}", url: "${env.JENKINS_JOB_URL}/job/$currentRelease/job/ConfirmReleaseConfigs/lastBuild/api/json"
					if (configStatusUrl.getStatus() == 404) {
						error("ConfirmReleaseConfigs has not yet been run. Please complete a successful build.")
					} else {
						def configStatusJson = new JsonSlurper().parseText(configStatusUrl.getContent())
						if (configStatusJson['result'] != "SUCCESS"){
							error("Most recent ConfirmReleaseConfigs build status: " + configStatusJson['result'] + ". Please complete a successful build.")
						}
					}
				}
			}
		}
		stage('Setup: Download Ortholog files from PANTHER'){
			steps{
				script{
					sh "wget ftp://ftp.pantherdb.org/ortholog/current_release/Orthologs_HCOP.tar.gz"
					sh "tar -xvf Orthologs_HCOP.tar.gz"
					sh "wget ftp://ftp.pantherdb.org/ortholog/current_release/QfO_Genome_Orthologs.tar.gz"
					sh "tar -xvf QfO_Genome_Orthologs.tar.gz"
				}
			}
		}
		stage('Setup: Download alternate ID mapping files'){
			steps{
				script{
					sh "wget -O mmus_alternate_ids.txt http://www.informatics.jax.org/downloads/reports/HGNC_homologene.rpt"
					sh "wget -O rnor_alternate_ids.txt ftp://ftp.rgd.mcw.edu/pub/data_release/GENES_RAT.txt"
					sh "wget -O xtro_alternate_ids.txt ftp://ftp.xenbase.org/pub/GenePageReports/GenePageEnsemblModelMapping.txt"
					sh "wget -O drer_alternate_ids.txt https://zfin.org/downloads/ensembl_1_to_1.txt"
				}
			}
		}			
		// This stage builds the jar file using maven.
		stage('Setup: Build jar file') {
			steps {
				script {
					sh "mvn clean compile assembly:single"
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
		// Logs and data files generated by this step are archived in the Reactome S3 bucket.
		// All files are then deleted on server.
		stage('Post: Archive Outputs'){
			steps{
				script{
					def s3Path = "${env.S3_RELEASE_DIRECTORY_URL}/${currentRelease}/orthopairs"
					sh "mkdir -p data/orthopairs/"
					sh "mv *alternate_ids.txt data/"
					sh "mv ${currentRelease}/* data/orthopairs/"
					sh "mv *.gz data/"
					sh "gzip -r logs/ data/"
					sh "aws s3 --no-progress --recursive cp logs/ $s3Path/logs/"
					sh "aws s3 --no-progress --recursive cp data/ $s3Path/data/"
					sh "rm -r ${currentRelease} logs data"
					sh "rm *Orthologs*"
				}
			}
		}
	}
}
