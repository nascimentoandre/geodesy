import numpy as np


def dms2degrees(coord):
    """
    Converts degrees, minutes and seconds to decimal degrees

    Parameters
    -----------
    coord: list or tuple
        Array containing a coordinate in degrees, minutes and seconds: [degrees, minutes, seconds]

    Returns
    ---------
    float
        Coordinate in decimal degrees
    """
    if coord[0] < 0:
        degrees = coord[0] - coord[1]/60 - coord[2]/3600
    else:
        degrees = coord[0] + coord[1]/60 + coord[2]/3600
    return degrees


def degrees2dms(coord):
    """
    Converts decimal degrees to degrees, minutes and seconds

    Parameters
    -----------
    coord: float
        Coordinate in decimal degrees

    Returns
    ---------
    list
        Degrees, minutes and seconds coordinate
    """
    coordstr = str(coord).split('.')[1]
    rest1 = float('0.'+coordstr)
    m = rest1 * 60
    minutes = float(str(m).split('.')[0])
    rest2 = float('0.'+str(m).split('.')[1])
    seconds = rest2 * 60
    return [float(str(coord).split('.')[0]), minutes, seconds]


def geod2cart(lamb, phi, h, elip):
    """
    Converts geodesic coordinates to three dimensional Cartesian coordinates

    Parameters
    -----------
    lamb: float
        Longitude in decimal degrees
    phi: float
        Latitude in decimal degrees
    h: float
        Geometric altitude in meters
    elip: object
        An instance of the class Ellipsoid

    Returns
    ---------
    float
        X component of the Cartesian coordinate
    float
        Y component of the Cartesian coordinate
    float
        Z component of the Cartesian coordinate
    """
    X = (elip.grandeNormal(phi) + h) * \
        np.cos(np.radians(phi)) * np.cos(np.radians(lamb))
    Y = (elip.grandeNormal(phi) + h) * \
        np.cos(np.radians(phi)) * np.sin(np.radians(lamb))
    Z = (elip.grandeNormal(phi)*(1 - elip.e1**2) + h) * np.sin(np.radians(phi))
    return X, Y, Z


def cart2geod(X, Y, Z, elip, dms=False):
    """
    Converts three dimensional Cartesian coordinates to geodesic coordinates

    Parameters
    -----------
    X: float
        X component of the Cartesian coordinate (in meters)
    Y: float
        Y component of the Cartesian coordinate (in meters)
    Z: float
        Z component of the Cartesian coordinate (in meters)
    elip: object
        An instance of the class Ellipsoid

    Returns
    ---------
    float
        Longitude
    float
        Latitude
    float
        Geometric Altitude in meters
    """
    u = np.arctan((Z/(np.sqrt(X**2+Y**2)))*(elip.a/elip.b))
    lamb = np.rad2deg(np.arctan(Y/X))
    phi = np.rad2deg(np.arctan((Z + elip.e2**2*elip.b*np.sin(u)**3)/((np.sqrt(X**2+Y**2)
                                                                      - elip.e1**2*elip.a*np.cos(u)**3))))
    h = (np.sqrt(X**2+Y**2)/np.cos(np.deg2rad(phi))) - elip.grandeNormal(phi)
    # Fixing lamb value
    if X < 0 and Y < 0:
        lamb = -(180-lamb)
    elif X < 0 and Y >= 0:
        lamb = 180 + lamb
    return lamb, phi, h
