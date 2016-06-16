#!/usr/bin/env python2

from math import sin, cos, radians, sqrt
import numpy as np
from os import makedirs, path
from subprocess import call, Popen, PIPE

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
    # MW, FLEN, DLEN, FWID, DWID, DTOP, STK, DIP, RAK, ELAT, ELON, SHYPO, DHYPO
    return lat, lon, depth, \
            Mw, fault_length, 0.1, fault_width, 0.1, depth_tcl, strike, dip, rake, \
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



if __name__ == "__main__":
    lats, lons, depths, MAG, FLEN, DLEN, FWID, DWID, DTOP, STK, DIP, RAK, ELAT, ELON, SHYPO, DHYPO = \
            CreateSimplifiedFiniteFaultFromFocalMechanism()

    GSF_DIR = 'Gsf'
    SRF_DIR = 'Srf'
    STOCH_DIR = 'Stoch'
    for folder in [GSF_DIR, SRF_DIR, STOCH_DIR]:
        if not path.exists(folder):
            makedirs(folder)

    GSF_BIN = '/hpc/home/rwg43/Bin/fault_seg2gsf'
    SRF_BIN = '/hpc/home/rwg43/Bin/genslip-v3.3'
    STOCH_BIN = '/hpc/home/rwg43/Bin/srf2stoch'
    VELFILE = 'lp_generic1d-gp01.vmod'

    NX = '%.0f' % (FLEN / DLEN)
    NY = '%.0f' % (FWID / DWID)

    # GSF 2 SRF
    DT = 0.025
    SEED = 1129571

    FILE_ROOT = 'm%.2f-%.1fx%.1f_s%d' % (MAG, FLEN, FWID, SEED)
    GSF_FILE = 'm%.2f-%.2fx%.2f.gsf' % (MAG, DLEN, DWID)
    SRF_FILE = '%s.srf' % (FILE_ROOT)
    STOCH_FILE = '%s.stoch' % (FILE_ROOT)

    # SRF 2 STOCH
    DX = 2.0
    DY = 2.0

    with open('%s/%s' % (GSF_DIR, GSF_FILE), 'w') as gsfp:
        gexec = Popen([GSF_BIN, 'read_slip_vals=0'], stdin = PIPE, stdout = gsfp)
        gexec.communicate('1\n%f %f %f %d %d %d %f %f %s %s' % \
                (ELON, ELAT, DTOP, STK, DIP, RAK, FLEN, FWID, NX, NY))

    with open('%s/%s' % (SRF_DIR, SRF_FILE), 'w') as srfp:
        call([SRF_BIN, 'read_erf=0', 'write_srf=1', 'read_gsf=1', 'write_gsf=0', \
                'infile=%s/%s' % (GSF_DIR, GSF_FILE), 'mag=%f' % (MAG), 'nx=%s' % (NX), \
                'ny=%s' % (NY), 'ns=1', 'nh=1', 'seed=%d' % (SEED), \
                'velfile=%s' % (VELFILE), 'shypo=%f' % (SHYPO), 'dhypo=%f' % (DHYPO), \
                'dt=%f' % DT, 'plane_header=1'], \
                stdout = srfp)

    with open('%s/%s' % (STOCH_DIR, STOCH_FILE), 'w') as stochp:
        with open('%s/%s' % (SRF_DIR, SRF_FILE), 'r') as srfp:
            call([STOCH_BIN, 'dx=%f' % (DX), 'dy=%f' % (DY)], stdin = srfp, stdout = stochp)

