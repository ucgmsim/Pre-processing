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
		export TMP=/tmp/`date|md5sum|cut -c-32`
		echo $TMP
		mkdir -p $TMP
		cd $TMP
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
		rm -rf /tmp/${env.HUDSON_SERVER_COOKIE}
		""" 
            }
        }
    }
}
