source /var/lib/jenkins/py3env/bin/activate
pip install -r requirements.txt

mkdir -p /tmp/${id}
cd /tmp/${id}
rm -rf qcore
git clone https://github.com/ucgmsim/qcore.git
pip install --no-deps ./qcore/
cd -
pytest --black --ignore=geoNet --ignore=NonUniformGrid --ignore=RegionalSeismicityTectonics --ignore=SrfGen/NHM/deprecated --ignore=test
if [ ! -z ${id} ]; then rm -rf /tmp/${id};echo '/tmp/${id} deleted';else echo "${id} is not set";fi;
