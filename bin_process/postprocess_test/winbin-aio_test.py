#!/usr/bin/env python2
"""
Tests winbin-aio bo comparing known outputs (in TestCaseVel) to ones generated.
Additional source file is passed to modify variables for testing outputs/inputs.
If the ASCII files aren't identicle:
    Compare headers (display if different) make sure it contains number of values.
    Compare vel readings, print number different, average and max relative differences.

@author Viktor Polak
@contact viktor.polak@canterbury.ac.nz
@date 3 May 2016

USAGE: execute from current directory only (winbin-aio.py must be in parent directory).
    '$ ./winbin-aio_test.py' or '$ python winbin-aio_test.py'

ISSUES: there should only be one main params file in parent, others should extend!
"""

from os.path import basename, exists
from subprocess import call
import hashlib
from decimal import *
from glob import glob

# smaller files only, read at once
def sha512sum(filepath):
    hasher = hashlib.sha512()
    with open(filepath, 'rb') as fp:
        hasher.update(fp.read())
    return hasher.hexdigest()



print('Running winbin-aio on local test case...')
with open('/dev/null', 'w') as sink:
    call(['python', 'winbin-aio.py', 'test_mode'], cwd='..', stderr=sink)

testcases = map(basename, glob('TestCaseVel/*'))
for case in testcases:
    if not exists('Vel/' + case):
        print('WARNING: ASCII file not created: ' + case)
        continue

    # hash check
    if sha512sum('TestCaseVel/' + case) == sha512sum('Vel/' + case):
        print('Success: Files are exact copies: ' + case)
        continue
    else:
        print('INVESTIGATING DIFFERENCES: ' + case)


    t_fp = open('TestCaseVel/' + case, 'r')
    r_fp = open('Vel/' + case, 'r')
    test_next = t_fp.readline
    result_next = r_fp.readline
    test_line = ''
    result_line = ''
    def read_next():
        global test_line, result_line
        test_line = test_next()
        result_line = result_next()
    def abort():
        global t_fp, r_fp
        t_fp.close()
        r_fp.close()

    # === header check (2 lines)
    for header_line in xrange(2):
        read_next()
        if test_line != result_line:
            print('header (line %d) diff:' %(header_line + 1))
            print('< ' + test_line)
            print('> ' + result_line)

    # number of values recorded
    try:
        n_test = int(test_line.split()[0])
        n_result = int(result_line.split()[0])
    except ValueError:
        print('Invalid header, number of values not found.')
        abort()
        continue
    if n_test != n_result:
        print('TestCaseVel has %d values but ResultVel has %d.' %(n_test, n_result))
        abort()
        continue

    # === value check
    read_next()
    # lines can have different no. of values, track which index in line we're at
    t_i = 0
    r_i = 0
    t_vals = test_line.split()
    r_vals = result_line.split()
    # determine significant places
    sig_fig = len(t_vals[0].split('e')[0].lstrip('+-').replace('.', ''))
    getcontext().prec = sig_fig

    def next_vals():
        global t_i, r_i, t_vals, r_vals, t_fp, r_fp
        t_i += 1
        r_i += 1
        if len(t_vals) == t_i:
            t_i = 0
            t_vals = t_fp.readline().split()
        if len(r_vals) == r_i:
            r_i = 0
            r_vals = r_fp.readline().split()

    n_diff = 0
    sum_rel = Decimal()
    rel_max = Decimal()
    rel_t = ''
    rel_r = ''
    for _ in xrange(n_test):
        t_vel = Decimal(t_vals[t_i])
        r_vel = Decimal(r_vals[r_i])
        if t_vel != r_vel:
            n_diff += 1
            diff = abs(t_vel - r_vel)
            rel = abs(diff / t_vel)
            sum_rel += rel
            if rel > rel_max:
                rel_max = rel
                rel_t = t_vals[t_i]
                rel_r = r_vals[r_i]
        next_vals()
    print('Differing values: %d/%d' % (n_diff, n_test))
    print('Relative diff AVG: %f' % (sum_rel / n_diff))
    print('Relative diff MAX: %f (%s vs %s)' % (rel_max, rel_t, rel_r))

    # closes files
    abort()
