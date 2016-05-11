#!/usr/bin/env python2
"""
Tests gen_ts by comparing known outputs (in TestCaseTSFiles) to ones generated.
Additional source file is passed to modify variables for testing outputs/inputs.
Check using hashes if files are the same.

@author Viktor Polak
@contact viktor.polak@canterbury.ac.nz
@date 4 May 2016

USAGE: execute from current directory only (gen_ts.py must be in parent directory).
    '$ ./gen_ts_test.py' or '$ python gen_ts_test.py'

ISSUES: 
"""

from os.path import basename, exists
from subprocess import call
import hashlib
from glob import glob

# smaller files only, read at once
def sha512sum(filepath):
    hasher = hashlib.sha512()
    with open(filepath, 'rb') as fp:
        hasher.update(fp.read())
    return hasher.hexdigest()



print('Running gen_ts on local test case...')
with open('/dev/null', 'w') as sink:
    call(['python', 'gen_ts.py', 'test_mode'], cwd='..', stderr=sink)

testcases = map(basename, glob('TestCaseTSFiles/*'))
for case in testcases:
    if not exists('TSFiles/' + case):
        print('WARNING: binary file not created: ' + case)
        continue

    # hash check
    if sha512sum('TestCaseTSFiles/' + case) == sha512sum('TSFiles/' + case):
        print('Success: Files are exact copies: ' + case)
        continue
    else:
        print('WARNING: files are different: ' + case)

