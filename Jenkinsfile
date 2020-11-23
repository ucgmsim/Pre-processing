pipeline {
    agent any 
    stages {
        stage('Install dependencies') {
            steps {
                echo 'Install dependencies on Jenkins server (maybe unnecessary if test runs inside Docker)'
		
		sh """
		pwd
		env
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
		pip install -r requirements.txt
		cd ${env.WORKSPACE}
		echo ${env.WORKSPACE}
		mkdir -p /tmp/${env.HUDSON_COOKIE}

		cd /tmp/${env.HUDSON_COOKIE}
		git clone https://github.com/ucgmsim/qcore.git
		pip install --no-deps ./qcore/
		cd ${env.WORKSPACE}/test
		
		pytest --black --ignore=geoNet --ignore=NonUniformGrid --ignore=RegionalSeismicityTectonics --ignore=SrfGen/NHM/deprecated --ignore=test
		"""
            }
        }
        stage('Publish Artifacts') {
            steps {
                echo 'Save the assemblies generated from the compilation' 
            }
        }
    }
}
