#!/bin/bash
set -e
DOWNLOADS_DIR = "downloads"

DIR=$(dirname "$(readlink -f "$0")") # Directory of the script -- allows the script to invoked from anywhere
mkdir -p $DIR/$DOWNLOADS_DIR
cd $DIR/$DOWNLOADS_DIR

## Download files
echo "Downloading ftp://ftp.pantherdb.org/ortholog/current_release/Orthologs_HCOP.tar.gz..."
wget ftp://ftp.pantherdb.org/ortholog/current_release/Orthologs_HCOP.tar.gz
tar -xvf Orthologs_HCOP.tar.gz

echo "Downloading ftp://ftp.pantherdb.org/ortholog/current_release/QfO_Genome_Orthologs.tar.gz..."
wget ftp://ftp.pantherdb.org/ortholog/current_release/QfO_Genome_Orthologs.tar.gz
tar -xvf QfO_Genome_Orthologs.tar.gz

echo "Downloading http://www.informatics.jax.org/downloads/reports/HGNC_homologene.rpt..."
wget -O mmus_alternate_ids.txt http://www.informatics.jax.org/downloads/reports/HGNC_homologene.rpt

echo "Downloading https://download.rgd.mcw.edu/data_release/GENES.RAT.txt"
wget -O rnor_alternate_ids.txt https://download.rgd.mcw.edu/data_release/GENES.RAT.txt

echo "Downloading https://ftp.xenbase.org/pub/GenePageReports/GenePageEnsemblModelMapping.txt..."
wget -O xtro_alternate_ids.txt https://ftp.xenbase.org/pub/GenePageReports/GenePageEnsemblModelMapping.txt

echo "Downloading https://zfin.org/downloads/ensembl_1_to_1.txt..."
wget -O drer_alternate_ids.txt https://zfin.org/downloads/ensembl_1_to_1.txt

cd $DIR

## Update repo
git pull
## Create a new jar file with the orthopairs code
mvn clean compile assembly:single

## Ensures the correct jar file is obtained regardless of the orthopairs project version
orthopairs_jar_file=$(ls target/orthopairs-*-jar-with-dependencies.jar)

## Run Orthopairs file generating script
echo "java -jar $orthopairs_jar_file"
java -jar $orthopairs_jar_file

echo "Removing large files..."
rm *Orthologs*_HCOP*

echo "Finished Orthopairs"
