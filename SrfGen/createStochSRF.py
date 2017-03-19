#!/usr/bin/env python2

from random import uniform, randint
from time import time

from setSrfParams import *
from createSRF import CreateSRF_ff, CreateSRF_multi

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

def CreateSRF_ffdStoch():
    # used for differentiating multiple runs
    # most significant digits 3+ years appart
    # least significant digit 10 seconds appart
    run_id = str(time())[2:9]

    # file containing variability info
    metainfo = '%s_%s.txt' % (M_NAME, run_id)
    # overwrite file / list variables printed
    with open('Srf/%s' % (metainfo), 'w') as of:
        for column in ['filename', 'seed', 'mag', 'flen', 'fwid']:
            of.write('%s\t' % column)
        of.write('\n')

    # create each scenario
    for ns in xrange(N_SCENARIOS):
        # scenario number as string
        nss = str(ns).zfill(4)
        # increment seed if wanted
        seed = SEED + SEED_INC * ns // N_SEED_INC

        # randomise MAGnitude
        if V_MAG[0]:
            m_mag = round(MAG + \
                    uniform(-V_MAG[0], V_MAG[0]), 2)
        else:
            m_mag = MAG
        # randomise FaultLENgth
        if V_FLEN[0]:
            # maximum multiple of DLEN to shift
            maxd = (FLEN * V_FLEN[0]) // DLEN
            m_flen = FLEN + randint(-maxd, maxd) * DLEN
        else:
            m_flen = FLEN
        # randomise FaultWIDth
        if V_FWID[0]:
            # maximum multiple of DWID to shift
            maxd = (FWID * V_FWID[0]) // DWID
            m_fwid = FWID + randint(-maxd, maxd) * DWID
        else:
            m_fwid = M_FWID[case]

        output = '%s_%s_%.4d.srf' % (M_NAME, run_id, ns)

        # run createSRF with randomised parameters
        CreateSRF_ff(LAT, LON, m_mag, STK, \
                RAK, DIP, DT, PREFIX, seed, RVFRAC, \
                ROUGH, SLIP_COV, m_flen, \
                DLEN, m_fwid, DWID, DTOP, \
                SHYPO, DHYPO, outroot = output, \
                genslip = GENSLIP)

        # append stoch data to info file
        with open('Srf/%s' % (metainfo), 'a') as of:
            # filename
            of.write('%s\t' % output)
            # seed
            of.write('%i\t' % seed)
            # magnitude
            of.write('%s\t' % param_as_string(m_mag))
            # fault length
            of.write('%s\t' % param_as_string(m_flen))
            # fault width
            of.write('%s\t' % param_as_string(m_fwid))
            # end of line / end of this srf file
            of.write('\n')

def CreateSRF_multiStoch():
    # used for differentiating multiple runs
    # most significant digits 3+ years appart
    # least significant digit 10 seconds appart
    run_id = str(time())[2:9]

    # file containing variability info
    metainfo = '%s_%s.txt' % (M_NAME, run_id)
    # overwrite file / list variables printed
    with open('Srf/%s' % (metainfo), 'w') as of:
        for column in ['filename', 'seed', 'nseg', 'mag', 'flen', 'fwid']:
            of.write('%s\t' % column)
        of.write('\n')

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
                c_flen = []
                for i, x in enumerate(M_FLEN[case]):
                    # maximum multiple of DLEN to shift
                    maxd = (x * V_FLEN[case]) // M_DLEN[case][i]
                    c_flen.append(x + randint(-maxd, maxd) * M_DLEN[case][i])
                m_flen.append(c_flen)
            else:
                m_flen.append(M_FLEN[case])
            # randomise FaultWIDth
            if V_FWID[case]:
                # maximum multiple of DWID to shift
                maxd = (M_FWID[case][0] * V_FWID[case]) // M_DWID[case][0]
                m_fwid.append([M_FWID[case][0] + \
                        randint(-maxd, maxd) * M_DWID[case][0]] \
                                        * len(M_FWID[case]))
            else:
                m_fwid.append(M_FWID[case])

        output = '%s_%s_%.4d.srf' % (M_NAME, run_id, ns)

        # run createSRF with randomised parameters
        CreateSRF_multi(M_NSEG, M_SEG_DELAY, m_mag, M_MOM, \
                M_RVFAC_SEG, M_GWID, M_RUP_DELAY, m_flen, \
                M_DLEN, m_fwid, M_DWID, M_DTOP, M_STK, \
                M_RAK, M_DIP, M_ELON, M_ELAT, M_SHYPO, \
                M_DHYPO, DT, seed, RVFRAC, ROUGH, SLIP_COV, \
                M_NAME, CASES, output = output, genslip = GENSLIP)

        # append stoch data to info file
        with open('Srf/%s' % (metainfo), 'a') as of:
            # filename
            of.write('%s\t' % output)
            # seed
            of.write('%i\t' % seed)
            # segment distribution
            of.write('%s\t' % param_as_string(M_NSEG))
            # magnitude
            of.write('%s\t' % param_as_string(m_mag))
            # fault length
            of.write('%s\t' % param_as_string(m_flen))
            # fault width
            of.write('%s\t' % param_as_string(m_fwid))
            # end of line / end of this srf file
            of.write('\n')

if __name__ == '__main__':
    if TYPE == 3:
        CreateSRF_ffdStoch()
    elif TYPE == 4:
        CreateSRF_multiStoch()
