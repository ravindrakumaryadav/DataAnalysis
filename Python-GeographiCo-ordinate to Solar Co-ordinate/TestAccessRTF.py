import cartopy
from spacepy.coordinates import Coords
from spacepy.time import Ticktock


def geo2smLL(qGeo, MJD):
    np, nd = qGeo.shape  # Np,2 (lon/lat)

    sphG = np.zeros((np, 3))  # Spherical geo, r/lat/lon [Re,deg,deg]
    sphG[:, 0] = 1.0
    sphG[:, 1] = qGeo[:, 1]
    sphG[:, 2] = qGeo[:, 0]
    MJDs = np.zeros(np)
    MJDs[:] = MJD

    invec = Coords(sphG, 'GEO', 'sph', use_irbem=False)
    invec.ticks = Ticktock(MJDs, dtype="MJD")
    outvec = invec.convert('SM', 'sph')

    qSM = np.zeros((np, 2))

    qSM[:, 0] = outvec.data[:, 2]
    qSM[:, 1] = outvec.data[:, 1]

    qSM_Fix = FixLine(qSM)
    return qSM_Fix


def FixLine(qLL):
    Np, Nd = qLL.shape  # Lon/lat

    for n in range(Np - 1):
        diff = qLL[n + 1, 0] - qLL[n, 0]

        if (diff > 180.0):
            qLL[n + 1, 0] = qLL[n + 1, 0] - 360
        if (diff < -180):
            qLL[n + 1, 0] = qLL[n + 1, 0] + 360
    return qLL


def AddCoasts(Ax, MJD, toplate):
    CLs = cartopy.feature.COASTLINE

    for geom in CLs.geometries():
        np, nCL = qGeo.shape  # Np,2 (lon/lat)

        qGeo = np.array(geom.coords[:])
        qSM = geo2smLL(qGeo, MJD)

        n = len(qSM[:, 0])

        if (n > nCL):
            Ax.plot(qSM[:, 0], qSM[:, 1], transform=toplate, color='slategrey', zorder=4)
