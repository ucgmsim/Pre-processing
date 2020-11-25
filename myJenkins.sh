source /var/lib/jenkins/py3env/bin/activate
echo "============================ DEPENDENCY INSTALLED ============================"
pip install -r requirements.txt
mkdir -p /tmp/${id}
cd /tmp/${id}
rm -rf qcore
git clone https://github.com/ucgmsim/qcore.git
pip install --no-deps ./qcore/
cd -
echo `pwd`
echo "================================= TEST BEGINS ================================="
pytest --black --ignore=geoNet --ignore=NonUniformGrid --ignore=RegionalSeismicityTectonics --ignore=SrfGen/NHM/deprecated --ignore=test
echo "=================================== CLEAN UP =================================="
if [ ! -z ${id} ]; then rm -rf /tmp/${id};echo '/tmp/${id} deleted';else echo "${id} is not set";fi;
