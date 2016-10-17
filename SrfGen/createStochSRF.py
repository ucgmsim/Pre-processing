#!/usr/bin/env python2

from random import uniform
from time import time

from setSrfParams import *
from createSRF import CreateSRF_multi

def param_as_string(param):
    # 1st case: single value
    if type(param).__name__ != 'list':
        return '%s' % (param)
    # 2nd case: 1D array
    if type(param[0]).__name__ != 'list':
        return ' '.join(map(str, param))
    # 3rd case: 2D array
    return ' '.join([' '.join(map(str, param[x])) \
            for x in xrange(len(param))])

def CreateSRF_multiStoch():
    # used for differentiating multiple runs
    # most significant digits 3+ years appart
    # least significant digit 10 seconds appart
    run_id = str(time())[2:9]

    # create each scenario
    for ns in xrange(N_SCENARIOS):
        # scenario number as string
        nss = str(ns).zfill(4)
        # increment seed if wanted
        seed = SEED + SEED_INC * ns // N_SEED_INC

        # cases can be randomised individually
        m_mag = []
        m_flen = []
        m_fwid = []
        for case in xrange(len(CASES)):
            # randomise MAGnitude
            if V_MAG[case]:
                m_mag.append(round(M_MAG[case] + \
                        uniform(-V_MAG[case], V_MAG[case]), 2))
            else:
                m_mag.append(M_MAG[case])
            # randomise FaultLENgth
            if V_FLEN[case]:
                m_flen.append([round(x + uniform(-x * V_FLEN[case], \
                        x * V_FLEN[case]), 2) \
                        for x in M_FLEN[case]])
            else:
                m_flen.append(M_FLEN[case])
            # randomise FaultWIDth
            if V_FWID[case]:
                m_fwid.append([M_FWID[case][0] + round( \
                        uniform(-M_FWID[case][0] * V_FWID[case], \
                                M_FWID[case][0] * V_FWID[case]), 2)] \
                                        * len(M_FWID[case]))
            else:
                m_fwid.append(M_FWID[case])

        output = '%s_%s_%.4d.srf' % (M_NAME, run_id, ns)

        # run createSRF with randomised parameters
        CreateSRF_multi(M_NSEG, M_SEG_DELAY, m_mag, M_MOM, \
                M_RVFAC_SEG, M_GWID, M_RUP_DELAY, m_flen, \
                M_DLEN, m_fwid, M_DWID, M_DTOP, M_STK, \
                M_RAK, M_DIP, M_ELON, M_ELAT, M_SHYPO, \
                M_DHYPO, DT, seed, M_NAME, CASES, \
                output)

        # append stoch data to end of SRF file
        with open('Srf/%s' % (output), 'a') as of:
            of.write('QCSTOCH\n')
            for param, value in {'NSEG':M_NSEG, 'SEG_DELAY':M_SEG_DELAY, \
                    'MAG':m_mag, 'MOM':M_MOM, 'RVFAC_SEG': M_RVFAC_SEG, \
                    'GWID':M_GWID, 'RUP_DELAY':M_RUP_DELAY, 'FLEN':m_flen, \
                    'DLEN':M_DLEN, 'FWID':m_fwid, 'DWID':M_DWID, \
                    'DTOP':M_DTOP, 'STK':M_STK, 'RAK':M_RAK, 'DIP':M_DIP, \
                    'ELON':M_ELON, 'ELAT':M_ELAT, 'SHYPO':M_SHYPO, \
                    'DHYPO':M_DHYPO, 'DT':DT, 'SEED':seed, \
                    'CASES':CASES}.iteritems():
                of.write(' %s\n' % (param))
                of.write('  %s\n' % param_as_string(value))

if __name__ == '__main__':
    if TYPE == 4:
        CreateSRF_multiStoch()
