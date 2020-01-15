#!/bin/bash
set -e

DIR=$(dirname "$(readlink -f "$0")") # Directory of the script -- allows the script to invoked from anywhere
cd $DIR

## Update repo
git pull
## Create new jar file with orthopairs code
mvn clean compile assembly:single

## Ensures the correct jar file is obtained regardless of orthopairs project version
orthopairs_jar_file=$(ls target/orthopairs-*-jar-with-dependencies.jar)

## Run Orthopairs file generating script
echo "java -jar $orthopairs_jar_file"
java -jar $orthopairs_jar_file

echo "Finished Orthopairs"
