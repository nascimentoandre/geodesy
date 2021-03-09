import numpy as np

def dist_sphere(phi1, lamb1, phi2, lamb2, R):
    """
    Calculates the distance between two points on a sphere

    Parameters
    -----------
    phi1: float
        Latitude of point 1
    lamb1: float
        Latitude of point 1
    phi2: float
        Latitude of point 2
    lamb2: float
        Longitude of point 2
    R: float or int
        Earth radius

    Returns
    --------
    float
        Distance between the two points over a sphere of radius R
    """
    phi1, lamb1 = np.deg2rad(phi1), np.deg2rad(lamb1)
    phi2, lamb2 = np.deg2rad(phi2), np.deg2rad(lamb2)
    S = np.arccos(np.sin(phi1)*np.sin(phi2)+np.cos(phi1)
                  * np.cos(phi2)*np.cos(lamb2-lamb1))
    Se = R*S
    return Se
