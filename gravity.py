import numpy as np
from coordinate_conv import dms2degrees


def normal_gravity(phi, elip=SIRGAS2000, dms=False):
    """
    @parâmetros: phi - Latitude do ponto, em graus.
                 elip - Elipsoide referente aos dados. Parâmetro optativo, o padrão é o SIRGAS2000.
                 dms - Graus, minutos e segundos (bool). Parâmetro optativo, o padrão é False.
    @retorna: gama - Gravidade teórica de pontos com esta latitude de entrada, em m/s².
    """
    if dms:
        phi = dms2degrees(phi)
    phi = np.deg2rad(phi)
    sin2_phi = np.sin(phi)**2

    gama_a = 9.8321863685  # m/s²
    gama_b = 9.7803267715  # m/s²

    k = (elip.b*gama_a - elip.a*gama_b)/(elip.a*gama_b)
    gama = gama_b*((1+k*sin2_phi)/(np.sqrt(1-elip.e1**2*sin2_phi)))
    return gama


def free_air_correction(phi, H, elip=SIRGAS2000, dms=False):
    """
    @parâmetros: phi - Latitude do ponto, em graus.
                 H - Altitude ortométrica do ponto, em metros.
                 elip - Elipsoide referente aos dados. Parâmetro optativo, o padrão é o SIRGAS2000.
                 dms - Graus, minutos e segundos (bool). Parâmetro optativo, o padrão é False.
    @retorna: fac - Correção ar livre, em m/s².
    """
    gama = normal_gravity(phi)  # m/s²
    sin2_phi = np.sin(phi)**2
    ang_vel = 7292115.0E-11  # rad/s
    GM = 3986005.0E8  # m³/s²
    m = (ang_vel**2*elip.a**2*elip.b)/GM
    fac = (2*gama/elip.a)*H*(1+elip.f+m-2*elip.f *
                             sin2_phi)-((3*gama*H**2)/(elip.a**2))
    # faa = faa*1.0E5 # retirar o comentário caso queira a correçao em mGal
    return fac


def free_air_anomaly(gobs, normal_grav, fac):
    """
    @parâmetros: gobs - Valor de gravidade observada, em mGal.
                 H - Valor de gravidade teórica do elipsoide, em m/s².
                 fac - Valor da correção ar livre do ponto, em m/s².
    @retorna: faa - Anomalia ar livre, em mGal.
    """
    gobs = gobs/1.0E5  # convertendo para m/s²
    faa = gobs - normal_grav + fac
    faa = faa * 1.0E5  # convertendo para mGal
    return faa


def bouguer_correction(H, ro=2670):
    """
    @parâmetros: H - Altitude ortométrica do ponto, em metros.
    @retorna: boug_c - Correção Bouguer, em m/s²
    """
    boug_c = 2*np.pi*ro*6.67408E-11*H
    return boug_c


def bouguer_anomaly(gobs, normal_grav, fac, bc):
    """
    @parâmetros: gobs - Valor de gravidade observada, em mGal.
                 H - Valor de gravidade teórica do elipsoide, em m/s².
                 fac - Valor da correção ar livre do ponto, em m/s².
                 bc - Valor da correção Bouguer do ponto, em m/s².
    @retorna: boug_a - Anomalia Bouguer, em mGal.
    """
    gobs = gobs/1.0E5  # convertendo para m/s²
    boug_a = gobs - normal_grav + fac - bc
    boug_a = boug_a * 1.0E5  # convertendo para mGal
    return boug_a


def C(df, i):
    """
    @parâmetros: df - DataFrame com os dados.
                 i - Número da linha do DataFrame na qual o ponto está.
    @retorna: C_f - Número de geopotencial do ponto.
    """
    C_f = 0
    for j in range(i):
        dn = df['Dn(i-1->i) (m)'][j+1]
        gp = df['Gravidade (mgal)'][j+1] / 1.0E5
        C_tmp = gp * dn
        C_f += C_tmp
    return C_f
