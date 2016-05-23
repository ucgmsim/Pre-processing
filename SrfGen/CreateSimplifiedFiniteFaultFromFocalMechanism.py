#!/usr/bin/env python2

from math import sin, cos, radians
from numpy import matrix, zeros, array

def CreateSimplifiedFiniteFaultFromFocalMechanism( \
        Lat = -43.5029, Lon = 172.8284, Depth = 4.0, Mw = 5.0, \
        strike = 54.0, rake = 137.0, dip = 75.0, MwScalingRel = 'HanksBakun2002'):
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
    faultGeometry = MwScalingRelation(Mw, MwScalingRel)

    # determine the number of fault patches needed
    Lx = 1.0
    Ly = 1.0
    Nx = int(round(faultGeometry / Lx))
    Ny = int(round(faultGeometry / Ly))

    # use cartesian coordinate system to define the along strike and downdip
    # locations taking the center of the fault plane as (x,y)=(0,0)
    xPos = [(x + 0.5) * Lx - Nx * Lx / 2.0 for x in xrange(0, Nx, 1)]
    yPos = [Ny * Ly / 2.0 - (x + 0.5) * Ly for x in xrange(0, Ny, 1)]

    # now use a coordinate transformation to go from fault plane to North and
    # East cartesian plane (again with (0,0) ==(0,0)  )
    yPosSurfProj = yPos * cos(radians(dip))
    rotAngleAzimuth = radians((strike - 90))
    RotMatrix = matrix([[cos(rotAngleAzimuth), sin(rotAngleAzimuth)], \
            [-sin(rotAngleAzimuth), cos(rotAngleAzimuth)]])
    eastLocRelative = zeros(shape = [Ny, Nx])
    northLocRelative = zeros(shape = [Ny, Nx])
    depthLocRelative = zeros(shape = [Ny, Nx])
    for i in xrange(Nx):
        for j in xrange(Ny):
            A = RotMatrix * matrix([[xPos(i)], [yPosSurfProj(j)]])
            eastLocRelative[j, i] = array(A).tolist()[0][0]
            northLocRelative[j, i] = array(A).tolist()[1][0]
            depthLocRelative[j, i] = -yPos[j] * sin(radians(dip))

    # now use a coordinate transformation to go from the North East cartesian
    # plane to the spherical earth WGS84 coordinate system
    EarthRadius = 6371.0072
    LengthOneDegreeLatitude = radians(EarthRadius)
    lat = zeros(shape = [Nx, Ny])
    lon = zeros(shape = [Nx, Ny])
    depth = zeros(shape = [Nx, Ny])
    for i in xrange(Nx):
        for j in xrange(Ny):
            lat[i, j] = latDatum + northLocRelative[i, j] / LengthOneDegreeLatitude
            lon[i, j] = lonDatum + (eastLocRelative[i, j] / LengthOneDegreeLatitude) \
                    * 1 / (cos(radians(latDatum)))
            depth[i, j] = max([depthDatum + depthLocRelative[i, j], 0])

    # outputs for emod3d format _________
    # now determine the location of the topcenter of the fault plane (top edge,
    # center point along strike) - see note at top for V3 changes.  "tcl" below
    # means "topcenter location"
    # -position in local coord system
    xPos_tcl = 0
    yPos_tcl = Ny * Ly / 2.0
    # convert to NE system
    yPosSurfProj_tcl = yPos_tcl * cos(radians(dip))
    A = RotMatrix * matrix([[xPos_tcl], [yPosSurfProj_tcl]])
    eastLocRelative_tcl = array(A).tolist()[0][0]
    northLocRelative_tcl = array(A).tolist()[1][0]
    depthLocRelative_tcl = -yPos_tcl * sin(radians(dip))
    # convert to Lat, Lon, depth
    lat_tcl = latDatum + northLocRelative_tcl / LengthOneDegreeLatitude
    lon_tcl = lonDatum + (eastLocRelative_tcl / LengthOneDegreeLatitude) * 1 / (cos(radians(latDatum)))
    depth_tcl = max([depthDatum + depthLocRelative_tcl, 0])

    # output format will depend on how it is used in python
    # copy of lat, lon, depth as 1ist most likely redundant?
    # lat, lon, depth, latAsList, lonAsList, depthAsList,
    # MW, FLEN, DLEN, FWID, DWID, DTOP, STK, DIP, RAK, ELAT, ELON, SHYPO, DHYPO
    return lat, lon, depth, array(lat).reshape(-1,).tolist(), \
            array(lon).reshape(-1,).tolist(), array(depth).reshape(-1,).tolist(), \
            Mw, faultGeometry, 0.1, faultGeometry, 0.1, depth_tcl, strike, dip, rake, \
            lat_tcl, lon_tcl, 0.00, faultGeometry / 2.0

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

    # only Length and Width of equiv. square (both sqrt(A)) were ever used
    return sqrt(A)

