pipeline {
    agent any 
    stages {
        stage('Install dependencies') {
            steps {
                echo 'Install dependencies on Jenkins server (maybe unnecessary if test runs inside Docker)'
		
		sh """
		pwd
		env
		source /var/lib/jenkins/py3env/bin/activate
		cd ${env.WORKSPACE}
		pip install -r requirements.txt
		echo ${currentBuild}
		mkdir -p /tmp/${env.ghprbActualCommit}
		cd /tmp/${env.ghprbActualCommit}
		rm -rf qcore
		git clone https://github.com/ucgmsim/qcore.git
		pip install --no-deps ./qcore/
		"""
            }
        }
        stage('Run regression tests') {
            steps {
                echo 'Run pytest' 
//		sh """
//		(if docker run is used, remove -it option)
//		docker run  -v /var/lib/jenkins/workspace/Pre-processing:/home/root/qcore sungeunbae/qcore-ubuntu-minimal bash -c "cd /home/root/qcore/;python setup.py install; cd qcore/test; pytest -s;"
		sh """

		source /var/lib/jenkins/py3env/bin/activate
		cd ${env.WORKSPACE}
		pytest --black --ignore=geoNet --ignore=NonUniformGrid --ignore=RegionalSeismicityTectonics --ignore=SrfGen/NHM/deprecated --ignore=test
		"""
            }
        }
        stage('Teardown') {
            steps {
                echo 'Clean up'
		sh """
		if [ ! -z ${env.ghprbActualCommit} ]; then rm -rf /tmp/${env.ghprbActualCommit};echo '/tmp/${env.ghprbActualCommit} deleted';else echo "${env.ghprbActualCommit} is not set";fi;
		""" 
            }
        }
    }
}
