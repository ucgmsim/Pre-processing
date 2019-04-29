""" Command to run this test: 'python -m pytest -v -s test_gen_coords.py' 
If the test passed it will delete the files in the output folder. 
Otherwise it would not delete files. 
 
 Instructions: Sample1 folder contains a sample output taken from hypocentre. Its path is noted in the readme file. In that path you will find the 
 vm_params.yaml along with other 5 output files. Use them as the benchmark files.If you want another sample to be tested,
 create a similar folder structure like sample1 and store the relevant files there (e.g:sample2). While running the test change sample1 to sample2

Just to run : py.test -s (or) python -m pytest -s -v test_gen_coords.py
To know the code coverage : py.test --cov=test_gen_coords.py
To know the test coverage :python -m pytest --cov ../../gen_coords.py test_gen_coords.py
"""

from qcore import shared
from glob import glob
import os
import shutil
import getpass
from datetime import datetime
import errno


PATH_TO_SAMPLE_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), "sample1")
PATH_TO_SAMPLE_OUTDIR = os.path.join(PATH_TO_SAMPLE_DIR, "output")
PATH_TO_SAMPLE_INPUT_DIR = os.path.join(PATH_TO_SAMPLE_DIR, "input")
INPUT_FILENAME = "vm_params.yaml"
# print "PATH_TO_SAMPLE_OUTDIR: ",PATH_TO_SAMPLE_OUTDIR


PATH_UNDER_TEST = os.path.join(
    "/home/",
    getpass.getuser(),
    ("tmp_test_gen_coords_" + "".join(str(datetime.now()).split())).replace(".", "_"),
).replace(":", "_")
# print "PATH_TO_NEW_OUTDIR: ", PATH_UNDER_TEST
# PATH_FOR_PRG_TOBE_TESTED = os.path.join(os.path.abspath(os.path.dirname(__file__)), "../gen_coords.py")

PATH_FOR_PRG_TOBE_TESTED = os.path.join(
    os.path.dirname(__file__), "../../SrfGen/gen_coords.py"
)
SYMLINK_PATH = os.path.join(
    os.path.abspath(os.path.dirname(PATH_FOR_PRG_TOBE_TESTED)), INPUT_FILENAME
)
# print "symbolic: **** ", SYMLINK_PATH


def setup_module(scope="module"):
    """ create a symbolic link for vm_params.yaml"""
    print("---------setup_module------------")
    try:
        os.mkdir(PATH_UNDER_TEST)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    sample_path = os.path.join(PATH_TO_SAMPLE_INPUT_DIR, INPUT_FILENAME)
    #    os.symlink(sample_path,SYMLINK_PATH)
    os.symlink(sample_path, os.path.join(PATH_UNDER_TEST, INPUT_FILENAME))


def test_gencoords():
    """ test qcore/gen_coords.py """
    print("---------test_gencoords------------")
    shared.exe("python " + PATH_FOR_PRG_TOBE_TESTED + " " + PATH_UNDER_TEST)
    ref_files = glob(os.path.join(PATH_TO_SAMPLE_OUTDIR, "*.100"))
    for ref in ref_files:
        cmd = "diff {} {}".format(
            ref, os.path.join(PATH_UNDER_TEST, os.path.basename(ref))
        )
        out, err = shared.exe(cmd)
        assert len(out) == 0 and len(err) == 0

    shutil.rmtree(PATH_UNDER_TEST)


def teardown_module():
    """ delete the symbolic link for vm_params.yaml"""
    print("---------teardown_module------------")
    if os.path.isfile(SYMLINK_PATH):
        os.remove(SYMLINK_PATH)
