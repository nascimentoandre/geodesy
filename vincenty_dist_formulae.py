import numpy as np

def inverse_problem(phi1, lamb1, phi2, lamb2, elip):
    """
    Given the coordinates of two points and an ellipsoid, this function
    calculates the azimuths alpha1 and alpha2, and the ellipsoidal distance s.

    Parameters
    -----------
    phi1: float
        Latitude of point 1 in degrees
    lamb1: float
        Longitude of point 1 in degrees
    phi2: float
        Latitude of point 2 in degrees
    lamb2: float
        Longitude of point 2 in degrees
    elip: object
        Instance of the Ellipsoid class

    Returns
    --------
    float
        Ellipsoidal distance in meters
    float
        Azimuth from 1 to 2 in degrees
    float
        Azimuth from 2 to 1 in degrees
    """
    phi1, phi2 = np.deg2rad(phi1), np.deg2rad(phi2)

    L = np.deg2rad(lamb2-lamb1)
    lamb = L
    U1 = np.arctan((1 - elip.f)*np.tan(phi1))
    U2 = np.arctan((1 - elip.f)*np.tan(phi2))
    i = 0

    while True:
        i += 1
        sin_sigma = np.sqrt((np.cos(U2)*np.sin(lamb))**2 + (np.cos(U1)*np.sin(U2) -
                                                            np.sin(U1)*np.cos(U2)*np.cos(lamb))**2)
        cos_sigma = np.sin(U1)*np.sin(U2) + np.cos(U1) * \
            np.cos(U2)*np.cos(lamb)
        sigma = np.arctan2(sin_sigma, cos_sigma)
        sin_alpha = (np.cos(U1)*np.cos(U2)*np.sin(lamb))/(sin_sigma)
        cos2_alpha = 1 - sin_alpha**2
        cos_dsigm = cos_sigma - ((2*np.sin(U1)*np.sin(U2))/(cos2_alpha))
        C = (elip.f/16)*cos2_alpha*(4+elip.f*(4 - 3*cos2_alpha))
        lamb_prev = lamb
        lamb = L + (1 - C)*elip.f*sin_alpha*(sigma + C*sin_sigma*(cos_dsigm + C*cos_sigma *
                                                                  (-1+2*cos_dsigm**2)))

        dif = np.abs(lamb-lamb_prev)
        if dif < 1.0E-12 or i >= 2000:
            break

    u2 = cos2_alpha * ((elip.a**2-elip.b**2)/elip.b**2)
    A = 1 + (u2/16384)*(4096+u2*(-768+u2*(320-175*u2)))
    B = (u2/1024)*(256+u2*(-128+u2*(74-47*u2)))
    del_sigma = B*sin_sigma*(cos_dsigm+0.25*B*(cos_sigma*(-1+2*cos_dsigm**2)-(1/6)*B
                                               * cos_dsigm*(-3+4*sin_sigma**2)*(-3+4*cos_dsigm**2)))
    s = elip.b*A*(sigma-del_sigma)
    alpha1 = np.rad2deg(np.arctan2(np.cos(U2)*np.sin(lamb), np.cos(U1)*np.sin(U2)
                                   - np.sin(U1)*np.cos(U2)*np.cos(lamb)))
    alpha2 = np.rad2deg(np.arctan2(np.cos(U1)*np.sin(lamb), -np.sin(U1)*np.cos(U2)
                                   + np.cos(U1)*np.sin(U2)*np.cos(lamb)))

    # Normalizing azimuths to 0..360:

    if alpha1 < 0:
        alpha1 = alpha1 + 360
    if alpha2 < 0:
        alpha2 = alpha2 + 360

    return s, alpha1, alpha2

def problema_direto(phi1, lamb1, alpha1, s, elip):
    """
    Given the coordinate of a point, azimuth alpha1, ellipsoidal distance
    and an ellipsoid, returns end point (phi2, lamb2) and azimuth alpha2

    Parameters
    ----------
    phi1: float
        Latitude of point 1 in degrees
    lamb1: float
        Longitude of point 1 in degrees
    alpha1: float
        Azimuth from 1 to 2 in degrees
    s: float
        Ellipsoidal distance in meters
    elip: object
        Instance of the Ellipsoid class

    Returns
    --------
    float
        Latitude of point 2 in degrees
    float
        Longitude of point 2 in degrees
    float
        Azimuth from 2 to 1 in degrees

    """
    phi1, lamb1, alpha1 = np.deg2rad(
        phi1), np.deg2rad(lamb1), np.deg2rad(alpha1)
    U1 = np.arctan((1-elip.f)*np.tan(phi1))
    sigma1 = np.arctan2(np.tan(U1), np.cos(alpha1))
    sin_alpha = np.cos(U1)*np.sin(alpha1)
    cos2_alpha = 1 - sin_alpha**2
    u2 = cos2_alpha*((elip.a**2-elip.b**2)/elip.b**2)
    A = 1 + (u2/16384)*(4096+u2*(-768+u2*(320-175*u2)))
    B = (u2/1024)*(256+u2*(-128+u2*(74-47*u2)))
    sigma = s/(elip.b * A)
    i = 0

    while True:
        i += 1
        sig_aux = sigma
        dois_sigma_m = 2*sigma1 + sigma
        del_sigma = B*np.sin(sigma)*(np.cos(dois_sigma_m)+0.25*B*(np.cos(sigma)*(-1+2*(np.cos(dois_sigma_m))**2)-(1/6)*B
                                                                  * np.cos(dois_sigma_m)*(-3+4*(np.sin(sigma))**2) *
                                                                  (-3+4*(np.cos(dois_sigma_m))**2)))

        sigma = (s/(elip.b*A)) + del_sigma
        if sigma-sig_aux < 1.0E-12 or i >= 2000:
            break

    phi2 = np.rad2deg(np.arctan2(np.sin(U1)*np.cos(sigma)+np.cos(U1)*np.sin(sigma)*np.cos(alpha1),
                                 (1-elip.f)*np.sqrt(sin_alpha**2+(np.sin(U1)*np.sin(sigma)-np.cos(U1) *
                                                                  np.cos(sigma)*np.cos(alpha1))**2)))
    lamb = np.arctan2(np.sin(sigma)*np.sin(alpha1),
                      np.cos(U1)*np.cos(sigma)-np.sin(U1)*np.sin(sigma)*np.cos(alpha1))
    C = (elip.f/16)*cos2_alpha*(4+elip.f*(4-3*cos2_alpha))
    L = lamb - (1-C)*elip.f*sin_alpha*(sigma+C*np.sin(sigma) *
                                       (np.cos(dois_sigma_m)+C*np.cos(sigma)*(-1+2*(np.cos(dois_sigma_m))**2)))
    lamb2 = np.rad2deg(lamb1 + L)
    alpha2 = np.rad2deg(np.arctan2(sin_alpha, -np.sin(U1)*np.sin(sigma)
                                   + np.cos(U1)*np.cos(sigma)*np.cos(alpha1)))

    # Normalizing the azimuth to 0..360:
    if alpha2 < 0:
        alpha2 = alpha2 + 360

    return phi2, lamb2, alpha2
