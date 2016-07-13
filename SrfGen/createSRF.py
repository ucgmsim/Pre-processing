#!/usr/bin/env python2

from math import exp, log, sin, cos, radians, sqrt
import numpy as np
from os import makedirs, path
from subprocess import call, Popen, PIPE

GSF_DIR = 'Gsf'
SRF_DIR = 'Srf'
STOCH_DIR = 'Stoch'

GSF_BIN = '/hpc/home/rwg43/Bin/fault_seg2gsf'
FF_SRF_BIN = '/hpc/home/rwg43/Bin/genslip-v3.3'
PS_SRF_BIN = '/hpc/home/rwg43/Bin/generic_slip2srf'
STOCH_BIN = '/hpc/home/rwg43/Bin/srf2stoch'
VELFILE = 'lp_generic1d-gp01.vmod'

mag2mom = lambda mw : exp(1.5 * (mw + 10.7) * log(10.0))
get_nx = lambda FLEN, DLEN : '%.0f' % (FLEN / DLEN)
get_ny = lambda FWID, DWID : '%.0f' % (FWID / DWID)
get_fileroot = lambda MAG, FLEN, FWID, seed : \
        'm%.2f-%.1fx%.1f_s%d' % (MAG, FLEN, FWID, seed)
get_gsfname = lambda MAG, DLEN, DWID : \
        'm%.2f-%.2fx%.2f.gsf' % (MAG, DLEN, DWID)


def CreateSimplifiedFiniteFaultFromFocalMechanism( \
        Lat = -43.5029, Lon = 172.8284, Depth = 4.0, Mw = 5.8, \
        strike = 54, rake = 137, dip = 75, MwScalingRel = 'BerrymanEtAl2002'):
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
    xPos = np.array([(x + 0.5) * Lx - Nx * Lx / 2.0 for x in xrange(0, Nx, 1)])
    yPos = np.array([Ny * Ly / 2.0 - (x + 0.5) * Ly for x in xrange(0, Ny, 1)])

    # now use a coordinate transformation to go from fault plane to North and
    # East cartesian plane (again with (0,0) ==(0,0)  )
    yPosSurfProj = yPos * cos(radians(dip))
    rotAngleAzimuth = radians((strike - 90))
    RotMatrix = np.array([[cos(rotAngleAzimuth), sin(rotAngleAzimuth)], \
            [-sin(rotAngleAzimuth), cos(rotAngleAzimuth)]])
    eastLocRelative = np.zeros(shape = [Ny, Nx])
    northLocRelative = np.zeros(shape = [Ny, Nx])
    depthLocRelative = np.zeros(shape = [Ny, Nx])
    for i in xrange(Nx):
        for j in xrange(Ny):
            A = np.dot(RotMatrix, [[xPos[i]], [yPosSurfProj[j]]])
            eastLocRelative[j, i] = np.array(A).tolist()[0][0]
            northLocRelative[j, i] = np.array(A).tolist()[1][0]
            depthLocRelative[j, i] = -yPos[j] * sin(radians(dip))

    # now use a coordinate transformation to go from the North East cartesian
    # plane to the spherical earth WGS84 coordinate system
    EarthRadius = 6371.0072
    LengthOneDegreeLatitude = radians(EarthRadius)
    lat = np.zeros(shape = [Nx, Ny])
    lon = np.zeros(shape = [Nx, Ny])
    depth = np.zeros(shape = [Nx, Ny])
    for i in xrange(Nx):
        for j in xrange(Ny):
            lat[i, j] = Lat + northLocRelative[i, j] / LengthOneDegreeLatitude
            lon[i, j] = Lon + (eastLocRelative[i, j] / LengthOneDegreeLatitude) \
                    * 1 / (cos(radians(Lat)))
            depth[i, j] = max([Depth + depthLocRelative[i, j], 0])

    # outputs for emod3d format _________
    # now determine the location of the topcenter of the fault plane (top edge,
    # center point along strike) - see note at top for V3 changes.  "tcl" below
    # means "topcenter location"
    # -position in local coord system
    xPos_tcl = 0
    yPos_tcl = Ny * Ly / 2.0
    # convert to NE system
    yPosSurfProj_tcl = yPos_tcl * cos(radians(dip))
    A = np.dot(RotMatrix, [[xPos_tcl], [yPosSurfProj_tcl]])
    eastLocRelative_tcl = A[0][0]
    northLocRelative_tcl = A[1][0]
    depthLocRelative_tcl = -yPos_tcl * sin(radians(dip))
    # convert to Lat, Lon, depth
    lat_tcl = Lat + northLocRelative_tcl / LengthOneDegreeLatitude
    lon_tcl = Lon + (eastLocRelative_tcl / LengthOneDegreeLatitude) * 1 / (cos(radians(Lat)))
    depth_tcl = max([Depth + depthLocRelative_tcl, 0])

    # lat, lon, depth,
    # MW, FLEN, DLEN, FWID, DWID, DTOP, ELAT, ELON, SHYPO, DHYPO
    return lat, lon, depth, \
            Mw, fault_length, 0.1, fault_width, 0.1, depth_tcl, \
            lat_tcl, lon_tcl, 0.00, fault_width / 2.0

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
    if not path.exists(GSF_DIR):
        makedirs(GSF_DIR)
    with open('%s/%s' % (GSF_DIR, gsf_file), 'w') as gsfp:
        gexec = Popen([GSF_BIN, 'read_slip_vals=0'], stdin = PIPE, stdout = gsfp)
        gexec.communicate('1\n%f %f %f %d %d %d %f %f %s %s' % \
                (lon, lat, dtop, strike, dip, rake, flen, fwid, nx, ny))

def gen_srf(srf_file, gsf_file, mw, dt, nx, ny, seed, shypo, dhypo):
    if not path.exists(SRF_DIR):
        makedirs(SRF_DIR)
    with open('%s/%s' % (SRF_DIR, srf_file), 'w') as srfp:
        call([FF_SRF_BIN, 'read_erf=0', 'write_srf=1', 'read_gsf=1', 'write_gsf=0', \
                'infile=%s/%s' % (GSF_DIR, gsf_file), 'mag=%f' % (mw), 'nx=%s' % (nx), \
                'ny=%s' % (ny), 'ns=1', 'nh=1', 'seed=%d' % (seed), \
                'velfile=%s' % (VELFILE), 'shypo=%f' % (shypo), 'dhypo=%f' % (dhypo), \
                'dt=%f' % dt, 'plane_header=1'], \
                stdout = srfp)

def gen_stoch(stoch_file, srf_file, dx = 2.0, dy = 2.0):
    if not path.exists(STOCH_DIR):
        makedirs(STOCH_DIR)
    with open('%s/%s' % (STOCH_DIR, stoch_file), 'w') as stochp:
        with open('%s/%s' % (SRF_DIR, srf_file), 'r') as srfp:
            call([STOCH_BIN, 'dx=%f' % (dx), 'dy=%f' % (dy)], stdin = srfp, stdout = stochp)

def CreateSRF_ps2ps(lat = -43.5871, lon = 172.5761, depth = 5.461, mw = -1, mom = -1, \
        strike = 246, rake = 159, dip = 84, prefix = '', stoch = True):
    """
    Must specify either magnitude or moment (mw, mom).
    """
    mw = 4.8

    # Vs, density at hypocenter
    VS = 3.20
    RHO = 2.44

    # source time function, sampling step
    STYPE = 'cos'
    RISE_TIME = 0.5
    DT = 0.005

    # fixed parameters
    TARGET_AREA_KM = -1
    TARGET_SLIP_CM = -1
    NSTK = 1
    NDIP = 1
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

    NAME = ('m%f' % (mw)).replace('.', 'pt').rstrip('0')
    GSF_FILE = '%s%s.gsf' % (prefix, NAME)
    SRF_FILE = '%s%s.srf' % (prefix, NAME)
    STOCH_FILE = '%s%s.stoch' % (prefix, NAME)

    FLEN = NSTK * float(d_xy) # dx
    FWID = NDIP * float(d_xy) # dy

    ###
    # GENERATE GSF
    if not path.exists(GSF_DIR):
        makedirs(GSF_DIR)
    with open('%s/%s' % (GSF_DIR, GSF_FILE), 'w') as gsfp:
        gsfp.write('# nstk= %d ndip= %d\n' % (NSTK, NDIP))
        gsfp.write('# flen= %10.4f fwid= %10.4f\n' % (FLEN, FWID))
        gsfp.write('# LON  LAT  DEP(km)  SUB_DX  SUB_DY  LOC_STK  LOC_DIP  LOC_RAKE  SLIP(cm)  INIT_TIME  SEG_NO\n')
        gsfp.write('%d\n' % (NSTK * NDIP))
        for j in xrange(0, NDIP):
            for i in xrange(1, NSTK + 1):
                gsfp.write('%11.5f %11.5f %8.4f %8.4f %8.4f %6.1f %6.1f %6.1f %8.2f %8.3f %3d\n' \
                        % (lon, lat, depth, float(d_xy), float(d_xy), strike, dip, rake, float(slip), RUPT, 0))

    ###
    # GENERATE SRF
    if not path.exists(SRF_DIR):
        makedirs(SRF_DIR)
    call([PS_SRF_BIN, 'infile=%s/%s' % (GSF_DIR, GSF_FILE), 'outfile=%s/%s' % (SRF_DIR, SRF_FILE), \
            'outbin=0', 'stype=%s' % (STYPE), 'dt=%f' % (DT), 'plane_header=1', \
            'risetime=%f' % (RISE_TIME), \
            'risetimefac=1.0', 'risetimedep=0.0'])

    ###
    # CONVERT TO STOCH
    if stoch:
        if not path.exists(STOCH_DIR):
            makedirs(STOCH_DIR)
        gen_stoch(STOCH_FILE, SRF_FILE, dx = 2.0, dy = 2.0)

def CreateSRF_cmt2ff(lat = -43.5029, lon = 172.8284, depth = 4.0, mw = 5.8, \
        strike = 54, rake = 137, dip = 75, mw_scaling_rel = 'BerrymanEtAl2002', \
        dt = 0.025, seed = 1129571, prefix = '', stoch = True):

    lats, lons, depths, MAG, FLEN, DLEN, FWID, DWID, DTOP, ELAT, ELON, SHYPO, DHYPO = \
            CreateSimplifiedFiniteFaultFromFocalMechanism(Lat = LAT, Lon = LON, \
            Depth = DEPTH, Mw = MAG, strike = STK, rake = RAK, dip = DIP, dipMwScalingRel = MWSR)

    NX = get_nx(FLEN, DLEN)
    NY = get_ny(FWID, DWID)

    FILE_ROOT = get_fileroot(MAG, FLEN, FWID, seed)
    GSF_FILE = '%s%s' % (prefix, get_gsfname(MAG, DLEN, DWID))
    SRF_FILE = '%s%s.srf' % (prefix, FILE_ROOT)
    STOCH_FILE = '%s%s.stoch' % (prefix, FILE_ROOT)

    gen_gsf(GSF_FILE, ELON, ELAT, DTOP, STK, DIP, RAK, FLEN, FWID, NX, NY)
    gen_srf(SRF_FILE, GSF_FILE, MAG, dt, NX, NY, seed, SHYPO, DHYPO)
    if stoch:
        gen_stoch(STOCH_FILE, SRF_FILE, dx = 2.0, dy = 2.0)

def CreateSRF_ffd2ff(lat = -43.5452, lon = 172.6971, mw = 6.2, \
        flen = 16.0, dlen = 0.10, fwid = 9.0, dwid = 0.10, dtop = 0.63, \
        shypo = -2.0, dhypo = 6.0, dt = 0.025, seed = 1129571, \
        strike = 59, rake = 128, dip = 69, prefix = '', stoch = True):

    NX = get_nx(flen, dlen)
    NY = get_ny(fwid, dwid)
    FILE_ROOT = get_fileroot(mw, flen, fwid, seed)
    GSF_FILE = '%s%s' % (prefix, get_gsfname(mw, dlen, dwid))
    SRF_FILE = '%s%s.srf' % (prefix, FILE_ROOT)
    STOCH_FILE = '%s%s.stoch' % (prefix, FILE_ROOT)

    gen_gsf(GSF_FILE, lon, lat, dtop, strike, dip, rake, flen, fwid, NX, NY)
    gen_srf(SRF_FILE, GSF_FILE, mw, dt, NX, NY, seed, shypo, dhypo)
    if stoch:
        gen_stoch(STOCH_FILE, SRF_FILE, dx = 2.0, dy = 2.0)

if __name__ == "__main__":
    from setSrfParams import *
    if TYPE == 1:
        # point source to point source srf
        CreateSRF_ps2ps(lat = LAT, lon = LON, depth = DEPTH, mw = MAG, mom = MOM, \
                strike = STK, rake = RAK, dip = DIP, prefix = PREFIX, stoch = True)
    elif TYPE == 2:
        # point source to finite fault srf
        CreateSRF_cmt2ff(lat = LAT, lon = LON, depth = DEPTH, mw = MAG, \
        dt = DT, prefix = PREFIX, stoch = True)
    elif TYPE == 3:
        # finite fault descriptor to finite fault srf
        CreateSRF_ffd2ff(lat = LAT, lon = LON, mw = MAG, flen = FLEN, dlen = DLEN, \
        fwid = FWID, dwid = DWID, dtop = DTOP, strike = STK, rake = RAK, dip = DIP, \
        prefix = PREFIX, stoch = True)
    else:
        print('Bad type of SRF generation specified. Check parameter file.')

