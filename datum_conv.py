from coordinate_conv import geod2cart, cart2geod

parameters = {}
with open("conversion_parameters.dat") as infile:
    for n, line in enumerate(infile):
        if n != 0:
            conv, dx, dy, dz = line.split()[0], float(line.split()[1]), float(line.split()[2]), float(line.split()[3])
            parameters[conv] = {'dx': dx, 'dy': dy, 'dz': dz}

def conv_geod_datum(lamb, phi, h, elip1, elip2, dms=False):
    """
    Converts coordinates between different geodetic datums. Conversion parameters defined until now:
    1) SIRGAS2000 -> SAD69
    2) SAD69 -> SIRGAS2000
    3) SICAD -> SIRGAS2000
    4) SIRGAS2000 -> SICAD

    Parameters
    -----------
    lamb: float
        Longitude in degrees
    phi: float
        Latitude in degrees
    h: float
        Geometric altitude in degrees
    elip1: object
        Instance of the Ellipsoid class related to the origin coordinates
    elip: object
        Instance of the Ellipsoid class related to the converted coordinates

    Returns
    ---------
    float
        Converted longitude in degrees
    float
        Converted latitude in degrees
    h
        Converted geometric altitude in meters
    """
    X, Y, Z = geod2cart(lamb, phi, h, elip1)
    ####### FIX THIS ############
    modo = elip1.nome + '2' + elip2.nome
    ##########################
    X2, Y2, Z2 = X + parameters[modo]['dx'], Y + \
        parameters[modo]['dy'], Z + parameters[modo]['dz']
    return cart2geod(X2, Y2, Z2, elip2)
