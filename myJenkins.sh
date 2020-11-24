source /var/lib/jenkins/py3env/bin/activate
cd ${env.WORKSPACE}
pip install -r requirements.txt

mkdir -p /tmp/${env.EXECUTOR_NUMBER}
cd /tmp/${env.EXECUTOR_NUMBER}
rm -rf qcore
git clone https://github.com/ucgmsim/qcore.git
pip install --no-deps ./qcore/
cd ${env.WORKSPACE}
pytest --black --ignore=geoNet --ignore=NonUniformGrid --ignore=RegionalSeismicityTectonics --ignore=SrfGen/NHM/deprecated --ignore=test
if [ ! -z ${env.EXECUTOR_NUMBER} ]; then rm -rf /tmp/${env.EXECUTOR_NUMBER};echo '/tmp/${env.EXECUTOR_NUMBER} deleted';else echo "${env.EXECUTOR_NUMBER} is not set";fi;
