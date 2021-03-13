import numpy as np

def problema_inverso(phi1, lamb1, phi2, lamb2, elip):
    """
    @parâmetros: phi1 - latitude do ponto 1
                 lamb1 - longitude do ponto 1
                 phi2 - latitude do ponto 2
                 lamb2 - longitude do ponto 2
                 elip - elipsoide ao longo do qual a distância será calculada
                 dms - graus, minutos e segundos (bool)

    @retorna: s - distância no elipsoide
              alpha1 - azimute de 1 para 2
              alpha2 - azimute de 2 para 1
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

    # Normalizando os azimutes para 0..360:

    if alpha1 < 0:
        alpha1 = alpha1 + 360
    if alpha2 < 0:
        alpha2 = alpha2 + 360

    return s, alpha1, alpha2

def problema_direto(phi1, lamb1, alpha1, s, elip):
    """
    @parâmetros: phi1 - latitude do ponto 1
                 lamb1 - longitude do ponto1
                 alpha1 - azimute de 1 para 2
                 s - distância no elipsoide
                 elip - elipsoide ao longo do qual a distância será calculada
                 dms - graus, minutos e segundos (bool)

    @retorna: phi2 - latitude do ponto 2
              lamb2 - longitude do ponto 2
              alpha2 - azimute de 2 para 1
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

    # Normalizando o azimute para 0..360:
    if alpha2 < 0:
        alpha2 = alpha2 + 360

    return phi2, lamb2, alpha2
