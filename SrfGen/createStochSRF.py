#!/usr/bin/env python2

from math import exp, log
import os
from random import uniform, randint
from time import time

from setSrfParams import *
from createSRF import CreateSRF_ff, CreateSRF_multi

mag2mom = lambda mw : exp(1.5 * (mw + 10.7) * log(10.0))
mom2mag = lambda mom : (2 / 3. * log(mom) / log(10.0)) - 10.7

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

def randomise_rupdelay(delay, var, deps):
    """
    Randomise rupture delay by variability 'var'
    such that dependencies 'dep' are adjusted as needed.
    delay: input delay times
    var: absolute variability
    deps: segment which triggers or None if independent
    """
    for i in xrange(len(delay)):
        if var[i] == 0:
            # segments have no variability
            continue

        # calculate, apply diff
        diff = uniform(-var[i], var[i])
        delay[i] += diff
        # propagate diff forwards only
        f = [i]
        while len(f) > 0:
            target = f.pop(0)
            for j in xrange(len(deps)):
                if deps[j] == target:
                    delay[j] += diff
                    f.append(j)
    return delay

def CreateSRF_multiStoch():
    # used for differentiating multiple runs
    # most significant digits 3+ years appart
    # least significant digit 10 seconds appart
    run_id = str(time())[2:9]

    out_dir = os.path.dirname(PREFIX)
    if out_dir == '':
        out_dir = '.'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # file containing variability info
    metainfo = '%s_%s.txt' % (PREFIX.rstrip('_'), run_id)
    # overwrite file / list variables printed
    with open(metainfo, 'w') as of:
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
        m0_tot = mag2mom(MW_TOTAL)
        # randomise time delay
        m_rdelay = randomise_rupdelay(M_RUP_DELAY, V_RDELAY, D_RDELAY)
        for case in xrange(len(CASES)):
            # randomise MAGnitude
            if V_MAG[case]:
                m_mag.append(mag2mom(M_MAG[case]) / m0_tot \
                        + uniform(-V_MAG[case], V_MAG[case]))
            else:
                m_mag.append(mag2mom(M_MAG[case]) / m0_tot)
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
        # make sure moment ratio is not negative
        m_mag = [max(m_mag[i], 0.000000001) for i in xrange(len(m_mag))]
        # normalise moment ratios
        m_mag = [m_mag[i] / sum(m_mag) for i in xrange(len(m_mag))]
        # convert back to magnitudes
        m_mag = [mom2mag(m_mag[i] * m0_tot) for i in xrange(len(m_mag))]

        output = '%s_%s_%.4d' % (PREFIX.rstrip('_'), run_id, ns)

        # run createSRF with randomised parameters
        CreateSRF_multi(M_NSEG, M_SEG_DELAY, m_mag, M_MOM, \
                M_RVFAC_SEG, M_GWID, m_rdelay, m_flen, \
                M_DLEN, m_fwid, M_DWID, M_DTOP, M_STK, \
                M_RAK, M_DIP, M_ELON, M_ELAT, M_SHYPO, \
                M_DHYPO, DT, seed, RVFRAC, ROUGH, SLIP_COV, \
                output, CASES, genslip = GENSLIP)

        # append stoch data to info file
        with open(metainfo, 'a') as of:
            # filename
            of.write('%s.srf\t' % output)
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
