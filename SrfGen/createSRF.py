#!/usr/bin/env python2

from math import exp, log, sin, cos, radians, sqrt
import numpy as np
from os import makedirs, path
from subprocess import call, Popen, PIPE

GSF_DIR = 'Gsf'
SRF_DIR = 'Srf'
STOCH_DIR = 'Stoch'

GSF_BIN = '/nesi/projects/nesi00213/tools/fault_seg2gsf'
FF_SRF_BIN = '/nesi/projects/nesi00213/tools/genslip-v3.3'
PS_SRF_BIN = '/nesi/projects/nesi00213/tools/generic_slip2srf'
STOCH_BIN = '/nesi/projects/nesi00213/tools/srf2stoch'
VELFILE = 'lp_generic1d-gp01.vmod'

mag2mom = lambda mw : exp(1.5 * (mw + 10.7) * log(10.0))
get_nx = lambda FLEN, DLEN : '%.0f' % (FLEN / DLEN)
get_ny = lambda FWID, DWID : '%.0f' % (FWID / DWID)
get_fileroot = lambda MAG, FLEN, FWID, seed : \
        'm%.2f-%.1fx%.1f_s%d' % (MAG, FLEN, FWID, seed)
get_gsfname = lambda MAG, DLEN, DWID : \
        'm%.2f-%.2fx%.2f.gsf' % (MAG, DLEN, DWID)

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

def CreateSRF_ps(lat = -43.5871, lon = 172.5761, depth = 5.461, mw = -1, mom = -1, \
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

def CreateSRF_ff(lat, lon, mw, strike, rake, dip, dt, prefix, seed, \
        flen = None, dlen = None, fwid = None, dwid = None, dtop = None, \
        shypo = None, dhypo = None, stoch = True, depth = None, mwsr = None,
        corners = True, corners_file = 'cnrs.txt'):
    """
    Create a Finite Fault SRF.
    Calculates flen, dlen... if not supplied given depth and mwsr are keywords.
    """

    # only given point source parameters? calculate rest
    if flen == None:
        lats, lons, depths, flen, dlen, fwid, dwid, dtop, elat, elon, shypo, dhypo = \
                focal_mechanism_2_finite_fault(lat, lon, depth, \
                        mw, strike, rake, dip, MWSR)[:]

    # write corners file
    corners = get_corners(lat, lon, flen, fwid, dip, strike)
    hypocentre = get_hypocentre(lat, lon, flen, fwid, shypo, dhypo, strike, dip)
    write_corners(corners_file, hypocentre, corners)

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
        CreateSRF_ps(LAT, LON, DEPTH, MAG, MOM, STK, RAK, DIP, PREFIX, stoch = True)
    elif TYPE == 2:
        # point source to finite fault srf
        CreateSRF_ff(LAT, LON, MAG, STK, RAK, DIP, DT, PREFIX, SEED, \
                depth = DEPTH, mwsr = MWSR, stoch = True, \
                corners = True, corners_file = CORNERS)
    elif TYPE == 3:
        # finite fault descriptor to finite fault srf
        CreateSRF_ff(LAT, LON, MAG, STK, RAK, DIP, DT, PREFIX, SEED, \
                FLEN, DLEN, FWID, DWID, DTOP, SHYPO, DHYPO, stoch = True, \
                corners = True, corners_file = CORNERS)
    else:
        print('Bad type of SRF generation specified. Check parameter file.')

