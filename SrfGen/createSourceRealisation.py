#!/usr/bin/env python2

from math import exp, log
import os
from random import uniform, randint, normalvariate
from time import time
import yaml

import argparse
from qcore import simulation_structure

# Uncomment if ffdStoch or multiStoch are needed, otherwise will break when enabled and attempting to use psStoch
#from setSrfParams import *
from createSRF import CreateSRF_ff, CreateSRF_multi, CreateSRF_ps

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

        # randomise hypocentre position - SHYPO -flen/2 -> flen/2
        shypo = SHYPO
        if V_HYPO[0]:
            mn = max( - FLEN / 2, SHYPO - V_HYPO[0])
            mx = min( FLEN / 2, SHYPO + V_HYPO[0])
            shypo = uniform(mn, mx)
        # randomise hypocentre posinion - HYPO 0 -> fwid
        dhypo = DHYPO
        if V_HYPO[1]:
            mn = max(0, DHYPO - V_HYPO[1])
            mx = min(FWID, DHYPO + V_HYPO[1])
            dhypo[0][0] = uniform(mn, mx)
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

        output = '%s_%s_%.4d' % (M_NAME, run_id, ns)

        # run createSRF with randomised parameters
        CreateSRF_ff(LAT, LON, m_mag, STK, \
                RAK, DIP, DT, output, seed, rvfrac = RVFRAC, \
                rough = ROUGH, slip_cov = SLIP_COV, flen = m_flen, \
                dlen = DLEN, fwid = m_fwid, dwid = DWID, dtop = DTOP, \
                shypo = shypo, dhypo = dhypo, \
                genslip = GENSLIP, stoch = STOCH)

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


def create_ps_realisation(out_dir, fault_name, lat, lon, depth, mw_mean, mom, strike, rake, dip,
                          uncertainty_file, n_realisations=50, additional_options={}, dt=0.005, vs=3.20, rho=2.44, target_area_km=None,
                          target_slip_cm=None, stype='cos', rise_time=0.5, init_time=0.0, silent=False):
    """
    Creates SRF files using random variables.
    Nominal values are to be passed in, while the probability characteristics are to be found in a yaml file that
    declares the distribution and probability characteristics for that distribution.
    Any values that are not to be passed to the CreateSRF function are to be passed in as a dict of additional_options
    and will be saved as a json file. If the keys of the additional_options are present in the yaml file they will also
    be perturbed with each realisation and saved as such
    :param lat:
    :param lon:
    :param depth:
    :param mw_mean:
    :param mom:
    :param strike:
    :param rake:
    :param dip:
    :param additional_options: A dictionary containing
    :param n_realisations:
    :param dt:
    :param srf_path:
    :param stoch_path:
    :param vs:
    :param rho:
    :param target_area_km:
    :param target_slip_cm:
    :param stype:
    :param rise_time:
    :param init_time:
    :param silent:
    """

    # Read yaml settings file
    with open(uncertainty_file) as yaml_file:
        yaml_settings = yaml.load(yaml_file)

    # Generate standard options dictionary

    standard_options = {'depth': depth,
                        'mw': mw_mean,
                        'mom': mom,
                        'strike': strike,
                        'rake': rake,
                        'dip': dip,
                        'vs': vs,
                        'rho': rho,
                        'rise_time': rise_time}

    for setting in yaml_settings:
        key = yaml_settings[setting].keys()[0]
        if 'mean' in yaml_settings[setting][key]:
            additional_options.update({setting: yaml_settings[setting][key]['mean']})

    perturbed_standard_options = dict(standard_options)
    perturbed_additional_options = dict(additional_options)

    for ns in xrange(n_realisations):

        realisation_name = simulation_structure.get_realisation_name(fault_name, ns)
        realisation_srf_path = os.path.join(out_dir, simulation_structure.get_srf_location(realisation_name))
        realisation_stoch_path = os.path.join(out_dir, simulation_structure.get_stoch_location(realisation_name))

        realisation_stoch_path = '/'+os.path.join(*realisation_stoch_path.split('/')[:-1])

        for key in yaml_settings.keys():

            # Load the correct dictionaries to be read from/written to
            if key in standard_options.keys():
                dict_to_read_from = standard_options
                dict_to_update = perturbed_standard_options
            elif key in additional_options.keys():
                dict_to_read_from = additional_options
                dict_to_update = perturbed_additional_options
            else:
                continue

            # Do this after checking the key will be used
            random_type = yaml_settings[key].keys()[0]
            random_value = yaml_settings[key][random_type]

            # Apply the random variable
            if random_type == 'uniform':
                dict_to_update[key] = uniform(dict_to_read_from[key] - random_value['halfrange'],
                                              dict_to_read_from[key] + random_value['halfrange'])
            elif random_type == 'normal':
                dict_to_update[key] = normalvariate(dict_to_read_from[key], random_value['std_dev'])
            elif random_type == 'uniform_relative':
                dict_to_update[key] = uniform(dict_to_read_from[key]*(1-random_value['scalefactor']),
                                              dict_to_read_from[key]*(1+random_value['scalefactor']))

        # Save the extra args to a json file
        additional_args_fname = os.path.join(out_dir, simulation_structure.get_source_params_location(realisation_name))
        with open(additional_args_fname, 'w') as yamlf:
            yaml.dump(perturbed_additional_options, yamlf)

        #print("Making srf with standard dict:")
        #print(perturbed_standard_options)
        #print("and additional dict:")
        #print(perturbed_additional_options)

        CreateSRF_ps(lat, lon, perturbed_standard_options['depth'], perturbed_standard_options['mw'], perturbed_standard_options['mom'], perturbed_standard_options['strike'],
                     perturbed_standard_options['rake'], perturbed_standard_options['dip'], dt=dt, prefix=realisation_srf_path[:-4], stoch=realisation_stoch_path,
                     vs=perturbed_standard_options['vs'], rho=perturbed_standard_options['rho'], target_area_km=target_area_km,
                     target_slip_cm=target_slip_cm, stype=stype, rise_time=perturbed_standard_options['rise_time'],
                     init_time=init_time, silent=silent)


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
        for column in ['filename', 'seed', 'nseg', 'mag', 'flen', 'fwid', 'rupdelay']:
            of.write('%s\t' % column)
        of.write('\n')

    # create each scenario
    for ns in xrange(N_SCENARIOS):
        # scenario number as string
        nss = str(ns).zfill(4)
        # increment seed if wanted
        seed = SEED + SEED_INC * ns // N_SEED_INC

        # randomise time delay
        m_rdelay = randomise_rupdelay(M_RUP_DELAY, V_RDELAY, D_RDELAY)
        # randomise hypocentre position - SHYPO -flen/2 -> flen/2
        m_shypo = list(M_SHYPO)
        if V_HYPO[0]:
            shypo = M_SHYPO[0][0]
            flen = M_FLEN[0][0]
            mn = max( - flen / 2, shypo - V_HYPO[0])
            mx = min( flen / 2, shypo + V_HYPO[0])
            m_shypo[0][0] = uniform(mn, mx)
        # randomise hypocentre posinion - HYPO 0 -> fwid
        m_dhypo = list(M_DHYPO)
        if V_HYPO[1]:
            dhypo = M_DHYPO[0][0]
            fwid = M_FWID[0][0]
            mn = max(0, dhypo - V_HYPO[1])
            mx = min(fwid, dhypo + V_HYPO[1])
            m_dhypo[0][0] = uniform(mn, mx)
        # cases can be randomised individually
        m_mag = []
        m_flen = []
        m_fwid = []
        m0_tot = sum([mag2mom(M_MAG[i]) for i in xrange(len(CASES))])
        m0_target = mag2mom(MW_TOTAL)
        for case in xrange(len(CASES)):
            # randomise MAGnitude
            if V_MAG[case]:
                maxd = mag2mom(M_MAG[case]) * V_MAG[case]
                m_mag.append((mag2mom(M_MAG[case]) \
                        + uniform(-maxd, maxd)) / m0_tot)
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
        # normalise moment ratios
        m_mag = [m_mag[i] / sum(m_mag) for i in xrange(len(m_mag))]
        # convert back to magnitudes
        m_mag = [mom2mag(m_mag[i] * m0_target) for i in xrange(len(m_mag))]

        output = '%s_%s_%.4d' % (PREFIX.rstrip('_'), run_id, ns)

        # run createSRF with randomised parameters
        CreateSRF_multi(M_NSEG, M_SEG_DELAY, m_mag, M_MOM, \
                M_RVFAC_SEG, M_GWID, m_rdelay, m_flen, \
                M_DLEN, m_fwid, M_DWID, M_DTOP, M_STK, \
                M_RAK, M_DIP, M_ELON, M_ELAT, m_shypo, \
                m_dhypo, DT, seed, output, CASES, rvfrac = RVFRAC, \
                rough = ROUGH, slip_cov = SLIP_COV, genslip = GENSLIP, \
                stoch = STOCH)

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
            # fault width
            of.write('%s\t' % param_as_string(m_rdelay))
            # end of line / end of this srf file
            of.write('\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("type", help="Type of earthquake to create SRF for")
    parser.add_argument("lat", help="latitude")
    parser.add_argument("lon", help="longitude")
    parser.add_argument("depth", help="depth")
    parser.add_argument("mw_mean", help="Mean magnitude")
    parser.add_argument("mom", help="Moment of quake. Will use mw_mean if <0")
    parser.add_argument("strike", help="Strike to be simulated")
    parser.add_argument("rake", help="Rake to be simulated")
    parser.add_argument("dip", help="Dip to be simulated")
    parser.add_argument("n_realisations", help="The number of realisations to generate SRFs for")
    args = parser.parse_args()

    if args.type == 1:
        create_ps_realisation(args.lat, args.lon, args.depth, args.mw_mean, args.mom, args.strike, args.rake, args.dip,
                              args.n_realisations)
    elif args.type == 3:
        CreateSRF_ffdStoch()
    elif args.type == 4:
        CreateSRF_multiStoch()
