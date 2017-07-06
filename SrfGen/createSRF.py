#!/usr/bin/env python2

from math import exp, log, log10, sin, cos, radians, sqrt
import os
from shutil import copyfile
from subprocess import call, Popen, PIPE
import sys

import numpy as np
sys.path.append(os.path.abspath(os.curdir)) #if there is a local srf_config, use it to override the default one

from srf_config import *

script_dir = os.path.dirname(os.path.abspath(__file__))
VELFILE = VELFILE.replace('<REPO>', '%s/' % (script_dir))

def srf_join(in1, in2, out):
    """
    like srf_join.c but that wasn't working properly
    """
    with open(in1, 'r') as f1:
        s1 = f1.readlines()
    with open(in2, 'r') as f2:
        s2 = f2.readlines()
    o = open(out, 'w')

    # must be the same version
    assert(s1[0] == s2[0])
    o.write(s1[0])

    # planes - header
    p1 = int(s1[1].split()[1])
    p2 = int(s2[1].split()[1])
    o.write('PLANE %d\n' % (p1 + p2))
    q1 = 2 + p1 * 2
    q2 = 2 + p2 * 2
    # planes - content
    o.write(''.join(s1[2 : q1]))
    o.write(''.join(s2[2 : q2]))

    # points - header
    o.write('POINTS %d\n' % (int(s1[q1].split()[1]) + int(s2[q2].split()[1])))
    # points - rest of file
    o.write(''.join(s1[q1 + 1 : ]))
    o.write(''.join(s2[q2 + 1 : ]))

    o.close()

mag2mom = lambda mw : exp(1.5 * (mw + 10.7) * log(10.0))
mom2mag = lambda mom : (2 / 3. * log(mom) / log(10.0)) - 10.7
get_nx = lambda FLEN, DLEN : '%.0f' % (FLEN / DLEN)
get_ny = lambda FWID, DWID : '%.0f' % (FWID / DWID)
get_fileroot = lambda MAG, FLEN, FWID, seed : \
        'm%.2f-%.1fx%.1f_s%d' % (MAG, FLEN, FWID, seed)
get_gsfname = lambda MAG, DLEN, DWID : \
        'm%.2f-%.2fx%.2f.gsf' % (MAG, DLEN, DWID)
# Leonard 2014 Relations
def leonard(rake, A):
    # if dip slip else strike slip
    if round(rake % 360 / 90.) % 2:
        return 4.19 + log10(A)
    else:
        return 4.18 + log10(A)

# by earth radius
one_deg_lat = radians(6371.0072)

# comment placed in corners file (above/below hypocenter)
# make sure lines end with '\n', especially the last one.
corners_header = ('> header line here for specifics \n\
> This is the standard input file format where the \
hypocenter is first then for each \n\
>Hypocenter (reference??) \n', \
'> Below are the corners \
(first point repeated as fifth to close box \n')

def write_corners(filename, hypocentre, corners):
    """
    Write a corners text file (used to plot faults).
    """
    with open(filename, 'w') as cf:
        # lines beginning with '>' are ignored
        cf.write(corners_header[0])
        cf.write('%f %f\n' % (hypocentre))
        cf.write(corners_header[1])
        # 0 1 - draw in order to close box
        # 2 3
        for i in [0, 1, 3, 2, 0]:
            cf.write('%f %f\n' % tuple(corners[i]))

def get_corners(lat, lon, flen, fwid, dip, strike):
    """
    Return Corners of a fault. Indexes are as below.
    0 1
    2 3
    """
    # use cartesian coordinate system to define the along strike and downdip
    # locations taking the center of the fault plane as (x,y)=(0,0)
    xPos = np.array([-flen / 2.0, flen / 2.0])
    yPos = np.array([-fwid, 0.0])
    Nx = len(xPos)
    Ny = len(yPos)

    # now use a coordinate transformation to go from fault plane to North and
    # East cartesian plane (again with (0,0) ==(0,0)  )
    yPosSurfProj = yPos * cos(radians(dip))
    azimuth = radians(strike - 90)
    RotMatrix = np.array([[cos(azimuth), sin(azimuth)], \
            [-sin(azimuth), cos(azimuth)]])
    eastLocRelative, northLocRelative = \
            np.dot(RotMatrix, [np.tile(xPos, Ny), \
                    yPosSurfProj.repeat(Ny)]).reshape(2, Nx, Ny)

    # now use a coordinate transformation to go from the North East cartesian
    # plane to the spherical earth WGS84 coordinate system
    lats = lat + northLocRelative / one_deg_lat
    lons = lon + (eastLocRelative / one_deg_lat) * 1 / cos(radians(lat))

    # joined (lon, lat) pairs
    return np.dstack((lons.flat, lats.flat))[0]

def get_hypocentre(lat, lon, flen, fwid, shypo, dhypo, strike, dip):
    """
    Same logic as corners, for a single point.
    """
    y_pos_surf_proj_hyp = -dhypo * cos(radians(dip))

    azimuth = radians(strike - 90)
    rot_matrix = np.array([[cos(azimuth), sin(azimuth)], \
            [-sin(azimuth), cos(azimuth)]])
    A = np.dot(rot_matrix, [[shypo], [y_pos_surf_proj_hyp]])

    lat_hyp = lat + A[1][0] / one_deg_lat
    lon_hyp = lon + (A[0][0] / one_deg_lat) * 1.0 / (cos(radians(lat)))

    return lon_hyp, lat_hyp

def focal_mechanism_2_finite_fault(Lat, Lon, Depth, Mw, \
        strike, rake, dip, MwScalingRel):
    """
    Purpose: To create a finite fault geometry based on centroid moment tensor
    solution information and magnitude scaling relationships.
    The use of such finite fault geometry is for first order computation of
    source-to-site distances.

    revisions
    v2 - reoriented fault plane to be consistent with along strike direction
    v3 - added 'argout_emod3d' output and associated commands to output details
         which are input in the finite fault model generation for EMOD3D graves methodology

    assumptions:
    1) The centroid moment tensor represents the centroid of the fault plane,
    located half way along strike and downdip
    2) The fault area is computed from the median of a Moment scaling
    relationship and thus is assumed to be a 'typical' stress drop event for
    the Mw-relationship used.
    Should only be used for small events where the downdip width is smaller
    than the seismogenic depth (i.e. so the assumption of a square fault plane
    is reasonable); depth is reset to zero if above ground values determined
    3) srike in deg measured clockwise from N; rake in deg measured
    anti-clockwise from the strike direction (and in the range
    -180<lambda<180; dip in deg measured downward from horizontal

    Input variables:
    Lat    - the latitude of the focal mechanism (-ve below equator)
    Lon    - the lon of the focal mech (-180<lon<180)
    Depth  - the depth of the focal mech (in km)
    Mw     - the moment magnitude from the Mw tensor soln
    strike - the strike of the focal mech (only 1)
    rake   - rake (not used, but carried forward)
    dip    - the dip angle of the fault plane

    output variables:
    lat       - the latitude of the subfault
    lon       - the lon of the subfault
    depth     - depth of the subfault
    lonAsList - as a 1D array
    """

    # get the fault geometry (square edge length)
    fault_length, fault_width = MwScalingRelation(Mw, MwScalingRel)

    # determine the number of fault patches needed
    Lx = 1.0
    Ly = 1.0
    Nx = int(round(fault_length / Lx))
    Ny = int(round(fault_width / Ly))

    # use cartesian coordinate system to define the along strike and downdip
    # locations taking the center of the fault plane as (x,y)=(0,0)
    xPos = np.array([(x + Lx / 2.0) * Lx - Nx * Lx / 2.0 for x in xrange(Nx)])
    yPos = np.array([Ny * Ly / 2.0 - (x + Ly / 2.0) * Ly for x in xrange(Ny)])

    # now use a coordinate transformation to go from fault plane to North and
    # East cartesian plane (again with (0,0) ==(0,0)  )
    yPosSurfProj = yPos * cos(radians(dip))
    rotAngleAzimuth = radians(strike - 90)
    RotMatrix = np.array([[cos(rotAngleAzimuth), sin(rotAngleAzimuth)], \
            [-sin(rotAngleAzimuth), cos(rotAngleAzimuth)]])
    eastLocRelative, northLocRelative = \
            np.dot(RotMatrix, [np.tile(xPos, Ny), \
                    yPosSurfProj.repeat(Ny)]).reshape(2, Nx, Ny)
    depthLocRelative = (-yPos * sin(radians(dip))).repeat(Ny).reshape((Nx, Ny))

    # now use a coordinate transformation to go from the North East cartesian
    # plane to the spherical earth WGS84 coordinate system
    lats = Lat + northLocRelative / one_deg_lat
    lons = Lon + (eastLocRelative / one_deg_lat) * 1 / cos(radians(Lat))
    depths = np.maximum(Depth + depthLocRelative, 0)

    # determine topcenter of the fault plane (top edge, center point along strike)
    # see note at top for V3 changes.  "tcl" means "topcenter location"
    xPos_tcl = 0
    yPos_tcl = Ny * Ly / 2.0
    # convert to NE system
    yPosSurfProj_tcl = yPos_tcl * cos(radians(dip))
    A = np.dot(RotMatrix, [[xPos_tcl], [yPosSurfProj_tcl]])
    eastLocRelative_tcl = A[0][0]
    northLocRelative_tcl = A[1][0]
    depthLocRelative_tcl = -yPos_tcl * sin(radians(dip))
    # convert to Lat, Lon, depth
    lat_tcl = Lat + northLocRelative_tcl / one_deg_lat
    lon_tcl = Lon + (eastLocRelative_tcl / one_deg_lat) * 1 / cos(radians(Lat))
    depth_tcl = max([Depth + depthLocRelative_tcl, 0])

    # FLEN, DLEN, FWID, DWID, DTOP, ELAT, ELON, SHYPO, DHYPO
    DLEN = 0.1
    DWID = 0.1
    SHYPO = 0.00
    DHYPO = fault_width / 2.0
    return lats, lons, depths, \
            fault_length, DLEN, fault_width, DWID, depth_tcl, \
            lat_tcl, lon_tcl, SHYPO, DHYPO

###########################################################################
def MwScalingRelation(Mw, MwScalingRel):
    """
    Return the fault Area from the Mw and a Mw Scaling relation.
    """
    if MwScalingRel == 'HanksBakun2002':
        try:
            # only small magnitude case
            assert(Mw <= 6.71)
        except AssertionError as e:
            e.args += ('Cannot use HanksAndBakun2002 equation for Mw > 6.71',)
            raise
        A = 10  ** (Mw - 3.98)

    elif MwScalingRel == 'BerrymanEtAl2002':
        A = 10 ** (Mw - 4.18)
        # i.e. set L=W=sqrt(A) in Eqn 1 in [1] Stirling, MW, Gerstenberger, M, Litchfield, N, McVerry, GH, Smith, WD, Pettinga, JR, Barnes, P. 2007.
        # updated probabilistic seismic hazard assessment for the Canterbury region, GNS Science Consultancy Report 2007/232, ECan Report Number U06/6. 58pp.

    elif MwScalingRel == 'VillamorEtAl2001':
        A = 10 ** (0.75 * (Mw - 3.39))
        # i.e. Eqn 2 in [1] Stirling, MW, Gerstenberger, M, Litchfield, N, McVerry, GH, Smith, WD, Pettinga, JR, Barnes, P. 2007. 
        # updated probabilistic seismic hazard assessment for the Canterbury region, GNS Science Consultancy Report 2007/232, ECan Report Number U06/6. 58pp.

    else:
        print('Invalid MwScalingRel. Exiting.')
        exit()

    # length, width
    return sqrt(A), sqrt(A)


def gen_gsf(gsf_file, lon, lat, dtop, strike, dip, rake, flen, fwid, nx, ny):
    with open(gsf_file, 'w') as gsfp:
        gexec = Popen([GSF_BIN, 'read_slip_vals=0'], \
                stdin = PIPE, stdout = gsfp)
        gexec.communicate('1\n%f %f %f %d %d %d %f %f %s %s' \
                % (lon, lat, dtop, strike, dip, rake, flen, fwid, nx, ny))
        gexec.wait()

def gen_srf(srf_file, gsf_file, mw, dt, nx, ny, seed, shypo, dhypo, \
            genslip = '3.3', rvfrac = None, rough = None, slip_cov = None):
    with open(srf_file, 'w') as srfp:
        genslip_bin = '%s-v%s' % (FF_SRF_BIN, genslip)
        if int(genslip[0]) < 5:
            xstk = 'nx'
            ydip = 'ny'
        else:
            xstk = 'nstk'
            ydip = 'ndip'
        cmd = [genslip_bin, 'read_erf=0', 'write_srf=1', 'read_gsf=1', \
                'write_gsf=0', 'infile=%s' % (gsf_file), \
                'mag=%f' % (mw), '%s=%s' % (xstk, nx), \
                '%s=%s' % (ydip, ny), 'ns=1', 'nh=1', 'seed=%d' % (seed), \
                'velfile=%s' % (VELFILE), 'shypo=%f' % (shypo), \
                'dhypo=%f' % (dhypo), 'dt=%f' % dt, 'plane_header=1', \
                'srf_version=2.0']
        if rvfrac != None:
            cmd.append('rvfrac=%s' % (rvfrac))
        if rough != None:
            cmd.append('alpha_rough=%s' % (rough))
        if slip_cov != None:
            cmd.append('slip_sigma=%s' % (slip_cov))
        print('Creating SRF:\n%s' % ' '.join(cmd))
        call(cmd, stdout = srfp)

def gen_stoch(stoch_file, srf_file, dx = 2.0, dy = 2.0):
    out_dir = os.path.dirname(stoch_file)
    if out_dir != '' and not os.path.exists(out_dir):
        os.makedirs(out_dir)
    with open(stoch_file, 'w') as stochp:
        with open(srf_file, 'r') as srfp:
            call([STOCH_BIN, 'dx=%f' % (dx), 'dy=%f' % (dy)], \
                    stdin = srfp, stdout = stochp)

def CreateSRF_ps(lat, lon, depth, mw, mom, \
        strike, rake, dip, dt, prefix = 'source', stoch = None):
    """
    Must specify either magnitude or moment (mw, mom).
    """
    # Vs, density at hypocentre
    VS = 3.20
    RHO = 2.44

    # source time function, sampling step
    STYPE = 'cos'
    RISE_TIME = 0.5
    DT = 0.005

    # fixed parameters
    TARGET_AREA_KM = -1
    TARGET_SLIP_CM = -1
    RUPT = 0.0

    if mom <= 0:
        # magnitude assumed to be set, calculate moment
        mom = mag2mom(mw)

    if TARGET_AREA_KM > 0:
        dd = sqrt(TARGET_AREA_KM)
        ss = (mom * 1.0e-20) / np.product([TARGET_AREA_KM, VS, VS, RHO])
    elif TARGET_SLIP_CM > 0:
        dd = sqrt(mom * 1.0e-20) / np.product([TARGET_SLIP_CM, VS, VS, RHO])
        ss = TARGET_SLIP_CM
    else:
        aa = exp(2.0 * log(mom) / 3.0 - 14.7 * log(10.0))
        dd = sqrt(aa)
        ss = (mom * 1.0e-20) / np.product([aa, VS, VS, RHO])
    d_xy, slip = ('%.5e %.3f' % (dd, ss)).split()

    if prefix[-1] == '_':
        prefix = '%s%s' % (prefix, ('m%f' % (mw)).replace('.', 'pt'))
    gsf_file = '%s.gsf' % (prefix)
    srf_file = '%s.srf' % (prefix)
    out_dir = os.path.dirname(srf_file)
    if out_dir != '' and not os.path.exists(out_dir):
        os.makedirs(out_dir)

    flen = float(d_xy)
    fwid = float(d_xy)
    print('FLEN: %s FWID: %s' % (flen, fwid))

    ###
    # GENERATE GSF
    with open(gsf_file, 'w') as gsfp:
        gsfp.write('# nstk= 1 ndip= 1\n')
        gsfp.write('# flen= %10.4f fwid= %10.4f\n' % (flen, fwid))
        gsfp.write('# LON  LAT  DEP(km)  SUB_DX  SUB_DY  LOC_STK  LOC_DIP  LOC_RAKE  SLIP(cm)  INIT_TIME  SEG_NO\n')
        gsfp.write('1\n')
        gsfp.write('%11.5f %11.5f %8.4f %8.4f %8.4f %6.1f %6.1f %6.1f %8.2f %8.3f %3d\n' \
                % (lon, lat, depth, float(d_xy), float(d_xy), \
                strike, dip, rake, float(slip), RUPT, 0))

    ###
    # GENERATE SRF
    call([PS_SRF_BIN, 'infile=%s' % (gsf_file), \
            'outfile=%s' % (srf_file), 'outbin=0', \
            'stype=%s' % (STYPE), 'dt=%f' % (DT), 'plane_header=1', \
            'risetime=%f' % (RISE_TIME), \
            'risetimefac=1.0', 'risetimedep=0.0'])

    ###
    # CONVERT TO STOCH
    if stoch != None:
        stoch_file = '%s/%s.stoch' % (stoch, os.path.basename(prefix))
        gen_stoch(stoch_file, srf_file, dx = 2.0, dy = 2.0)

    # location of resulting SRF file
    return srf_file

def CreateSRF_ff(lat, lon, mw, strike, rake, dip, dt, prefix0, seed, \
        flen = None, dlen = None, fwid = None, dwid = None, \
        dtop = None, shypo = None, dhypo = None, stoch = None, depth = None, \
        mwsr = None, corners = True, corners_file = 'cnrs.txt', \
        genslip = '3.3', rvfrac = None, rough = None, slip_cov = None):
    """
    Create a Finite Fault SRF.
    Calculates flen, dlen... if not supplied given depth and mwsr are keywords.
    NOTE: depth is focal mech depth and only used when calculating flen etc...
    """

    # do not change input variables
    prefix = prefix0

    # only given point source parameters? calculate rest using scaling relation
    if flen == None:
        flen, dlen, fwid, dwid, dtop, elat, elon, shypo, dhypo = \
                focal_mechanism_2_finite_fault(lat, lon, depth, \
                        mw, strike, rake, dip, MWSR)[3:]

    # write corners file
    corners = get_corners(lat, lon, flen, fwid, dip, strike)
    hypocentre = get_hypocentre(lat, lon, flen, fwid, shypo, dhypo, strike, dip)
    write_corners(corners_file, hypocentre, corners)

    nx = get_nx(flen, dlen)
    ny = get_ny(fwid, dwid)
    if prefix[-1] == '_':
        prefix = '%s%s' % (prefix, get_fileroot(mw, flen, fwid, seed))
    gsf_file = '%s.gsf' % (prefix)
    srf_file = '%s.srf' % (prefix)
    out_dir = os.path.dirname(srf_file)
    if out_dir != '' and not os.path.exists(out_dir):
        os.makedirs(out_dir)

    gen_gsf(gsf_file, lon, lat, dtop, strike, dip, rake, flen, fwid, nx, ny)
    gen_srf(srf_file, gsf_file, mw, dt, nx, ny, seed, shypo, dhypo, \
            genslip = genslip, rvfrac = rvfrac, slip_cov = slip_cov, \
            rough = rough)
    if stoch != None:
        stoch_file = '%s/%s.stoch' % (stoch, os.path.basename(prefix))
        gen_stoch(stoch_file, srf_file, dx = 2.0, dy = 2.0)

    # print leonard Mw from A (SCR)
    print('Leonard 2014 Mw: %s' % (leonard(rake, fwid * flen)))
    # location of resulting SRF
    return srf_file

def CreateSRF_multi(nseg, seg_delay, mag0, mom0, rvfac_seg, gwid, rup_delay, \
        flen, dlen, fwid, dwid, dtop, stk, rak, dip, elon, elat, \
        shypo, dhypo, dt, seed, prefix0, cases, genslip = '3.3', \
        rvfrac = None, rough = None, slip_cov = None, stoch = None):

    # do not change any pointers
    mag = list(mag0)
    mom = list(mom0)
    prefix = prefix0

    genslip_bin = '%s-v%s' % (FF_SRF_BIN, genslip)
    if int(genslip[0]) < 5:
        xstk = 'nx'
        ydip = 'ny'
        rup_name = 'rup_delay'
    else:
        xstk = 'nstk'
        ydip = 'ndip'
        rup_name = 'rupture_delay'
    out_dir = os.path.dirname(prefix)
    if out_dir != '' and not os.path.exists(out_dir):
        os.makedirs(out_dir)

    casefiles = []
    for c, case in enumerate(cases):

        if mom[c] <= 0:
            mom[c] = mag2mom(mag[c])
        else:
            mag[c] = mom2mag(mom[c])

        nx = map(int, [get_nx(flen[c][f], dlen[c][f]) \
                for f in xrange(nseg[c])])
        ny = map(int, [get_ny(fwid[c][f], dwid[c][f]) \
                for f in xrange(nseg[c])])
        nx_tot = sum(nx)
        flen_tot = sum(flen[c])
        # hypocentre is given based on the first segment
        # needs to be relative to full length instead of first segment
        shyp_tot = shypo[c][0] + 0.5 * flen[c][0] - 0.5 * flen_tot

        xseg = '-1'
        if int(nseg[c] > 1):
            xseg = []
            sbound = 0.0
            for g in xrange(nseg[c]):
                sbound += flen[c][g]
                xseg.append(sbound - 0.5 * flen_tot)
            xseg = ','.join(map(str, xseg))

        case_root = '%s_%s' % (prefix.rstrip('_'), case)
        gsf_file = '%s.gsf' % (case_root)
        srf_file = '%s.srf' % (case_root)
        fsg_file = '%s_seg.in' % (case_root)
        casefiles.append(srf_file)

        with open(fsg_file, 'w') as fs:
            fs.write('%d\n' % (nseg[c]))
            for f in xrange(nseg[c]):
                fs.write('%f %f %f %.4f %.4f %.4f %.4f %.4f %i %i\n' % ( \
                        elon[c][f], elat[c][f], dtop[c][f], \
                        stk[c][f], dip[c][f], rak[c][f], \
                        flen[c][f], fwid[c][f], nx[f], ny[f]))
        # create GSF file
        call([GSF_BIN, 'read_slip_vals=0', 'infile=%s' % (fsg_file), \
                'outfile=%s' % (gsf_file)])
        os.remove(fsg_file)

        # calling with outfile parameter will:
        # append 's0000-h0000' to filename (v3.3)
        # ns = slip realisations, nh = hypocentre realisations
        with open(srf_file, 'w') as srfp:
            cmd = [genslip_bin, 'read_erf=0', 'write_srf=1', \
                    'seg_delay=%s' % (seg_delay[c]), 'read_gsf=1', \
                    'write_gsf=1', 'infile=%s' % (gsf_file), \
                    'nseg=%s' % (nseg[c]), 'nseg_bounds=%d' % (nseg[c] - 1), \
                    'xseg=%s' % (xseg), 'rvfac_seg=%s' % (rvfac_seg[c]), \
                    'gwid=%s' % (gwid[c]), 'mag=%f' % (mag[c]), \
                    '%s=%f' % (xstk, nx_tot), '%s=%f' % (ydip, ny[0]), \
                    'ns=1', 'nh=1', 'seed=%d' % (seed), \
                    'velfile=%s' % (VELFILE), 'shypo=%f' % (shyp_tot), \
                    'dhypo=%f' % (dhypo[c][0]), 'dt=%f' % (dt), \
                    'plane_header=1', 'side_taper=0.02', 'bot_taper=0.02', \
                    'top_taper=0.0', '%s=%s' % (rup_name, rup_delay[c]), \
                    'srf_version=2.0']
            if rvfrac != None:
                cmd.append('rvfrac=%s' % (rvfrac))
            if rough != None:
                cmd.append('alpha_rough=%s' % (rough))
            if slip_cov != None:
                cmd.append('slip_sigma=%s' % (slip_cov))
            call(cmd, stdout = srfp)
        #os.remove(gsf_file)
        # print leonard Mw from A (SCR)
        print('Leonard 2014 Mw: %s' % (leonard(rak[c][f], fwid[c][f] * flen[c][f])))

    # joined filename
    if prefix[-1] == '_':
        prefix = '%ss%d' % (prefix, seed)
    joined_srf = '%s.srf' % (prefix)
    # joined case files stored in separate file
    copyfile(casefiles[0], joined_srf)
    # join rest of cases into the first
    for i in xrange(len(cases) - 1):
        srf_join(joined_srf, casefiles[i + 1], joined_srf)
    # remove casefiles
    for casefile in casefiles:
        os.remove(casefile)
    if stoch != None:
        stoch_file = '%s/%s.stoch' % (stoch, os.path.basename(prefix))
        gen_stoch(stoch_file, joined_srf, dx = 2.0, dy = 2.0)

    # path to resulting SRF
    return joined_srf


if __name__ == "__main__":
    try:
        from setSrfParams import *
    except ImportError:
        print('setSrfParams.py not found.')
        exit(1)
    except:
        print "Error: setSrfParams.py has issues."
        print sys.exc_info()[0]
        raise

    # incomplete defaults compatibility
    try:
        RVFRAC
    except NameError:
        RVFRAC = None
    try:
        ROUGH
    except NameError:
        ROUGH = None
    try:
        SLIP_COV
    except NameError:
        SLIP_COV = None
    try:
        STOCH
    except NameError:
        STOCH = None


    if TYPE == 1:
        # point source to point source srf
        srf = CreateSRF_ps(LAT, LON, DEPTH, MAG, MOM, STK, RAK, DIP, DT, \
        PREFIX, stoch = STOCH)
    elif TYPE == 2:
        # point source to finite fault srf
        srf = CreateSRF_ff(LAT, LON, MAG, STK, RAK, DIP, DT, PREFIX, \
                SEED, depth = DEPTH, mwsr = MWSR, \
                stoch = STOCH, corners = True, \
                genslip = GENSLIP, rvfrac = RVFRAC, slip_cov = SLIP_COV, \
                rough = ROUGH)
    elif TYPE == 3:
        # finite fault descriptor to finite fault srf
        srf = CreateSRF_ff(LAT, LON, MAG, STK, RAK, DIP, DT, PREFIX, \
                SEED, FLEN, DLEN, FWID, DWID, DTOP, \
                SHYPO, DHYPO, stoch = STOCH, corners = True, \
                genslip = GENSLIP, rvfrac = RVFRAC, slip_cov = SLIP_COV, \
                rough = ROUGH)
    elif TYPE == 4:
        # multi segment finite fault srf
        srf = CreateSRF_multi(M_NSEG, M_SEG_DELAY, M_MAG, M_MOM, \
                M_RVFAC_SEG, M_GWID, M_RUP_DELAY, M_FLEN, \
                M_DLEN, M_FWID, M_DWID, M_DTOP, M_STK, \
                M_RAK, M_DIP, M_ELON, M_ELAT, M_SHYPO, \
                M_DHYPO, DT, SEED, PREFIX, CASES, \
                genslip = GENSLIP, rvfrac = RVFRAC, slip_cov = SLIP_COV, \
                stoch = STOCH, rough = ROUGH)
    else:
        print('Bad type of SRF generation specified. Check parameter file.')
        exit(1)

    # run SRF plots
    if len(sys.argv) > 1 and sys.argv[1] == 'noplot':
        exit(0)

    if TYPE == 1:
        print('Point source plotting not yet implemented.')
        exit(0)

    # start plotting
    from subprocess import call, Popen, PIPE

    for script in ['plot_srf_square.py', 'plot_srf_map.py']:
        path_tester = Popen(['which', script], stdout = PIPE)
        script_path = path_tester.communicate()[0].rstrip()
        path_tester.wait()
        if not os.path.exists(script_path):
            print('Plotting script %s is not available in PATH!' % (script))
            exit(1)

    print('Plotting SRF as square plot...')
    call(['plot_srf_square.py', srf])
    print('Plotting SRF as map plot...')
    call(['plot_srf_map.py', srf])
