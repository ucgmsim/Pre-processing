#!/usr/bin/env python2
"""
Converts something into something else.

@authors: Viktor Polak
@date 05 Aug 2016

ISSUES:
    could get HH from params.py but to reduce dependencies...
    could get it from a path/filename (depending how this is used)
"""

INFILE = 'lp_generic1d-gp01.vmod'
HH = '0.100'
OUTFILE = 'lp_generic1d-gp01.fd-h%s-pytest' % (HH)

# INFILE column index definitions
dd = 0
vp = 1
vs = 2
dn = 3
qp = 4
qs = 5

with open(INFILE, 'r') as inp:
    # first line = number of lines to follow
    in_items = int(inp.readline())
    # lines (1st dimention), values (2nd dimention)
    in_data = map(lambda l : map(float, l), \
            map(str.split, inp.readlines()))

k = lambda dd_sum : int(dd_sum / float(HH) + 1)

with open(OUTFILE, 'w') as outp:
    # write header
    outp.write('DEF HST\n')

    dd_sum = 0.0
    k_old = 1
    # reason for indexing instead of looping over lines:
    # TODO: use a 'k_next' to clean up write statements
    for i in xrange(in_items):
        # details in line
        dd_sum += in_data[i][dd]
        k_new = k(dd_sum)
        if k_new > k_old:
            k_old = k_new
            outp.write(output)

        # only stored if next k > k or last item
        output = ('%8.4f %8.4f %8.4f %10.2f %10.2f %9.3f\n' % ( \
                in_data[i][vp], in_data[i][vs], in_data[i][dn], \
                in_data[i][qp], in_data[i][qs], \
                # ZZ: (means what?)
                float(HH) * int(dd_sum / float(HH) + 0.5) + 0.5 * float(HH) \
                ))
    # loop doesn't check k_next
    outp.write(output)
