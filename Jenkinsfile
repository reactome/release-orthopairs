pipeline {
    agent any

    stages {
        stage('Build') {
            steps {
		script {
                    dir ('orthopairs') {
                  	sh 'mvn clean compile assembly:single'
                  	sh 'java -jar target/orthopairs-0.0.1-SNAPSHOT-jar-with-dependencies.jar'
                    }
          	}
            }
        }
    }
}
