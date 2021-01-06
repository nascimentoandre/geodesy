import numpy as np
from coordinate_conv import dms2degrees


def calc_comp_desvio(lamb_g, phi_g, lamb_a, phi_a, dms=False):
    """
    calc_comp_desvio(lamb_g, phi_g, lamb_a, phi_a, dms=False)

    Dadas as coordenadas geodésicas e astronômicas de um ponto, calcula as componentes meridiana e primeiro
    vertical do desvio vertical.

    Parâmetros
    ------------
    lamb_g: Longitude geodésica.
    phi_g: Latitude geodésica.
    lamb_a: Longitude astronômica.
    phi_a: Latitude astronômica.
    dms: Formato das coordenadas de entrada. Caso as coordenadas estejam em graus, minutos e segundos,
    deve-se atribuir um valor True.

    Retorna
    --------------
    comp_meridiana: Valor da componente meridiana em segundos.
    comp_primeiro_vertical: Valor da componente primeiro vertical em segundos.
    """
    if dms:
        lamb_g = dms2degrees(lamb_g)
        phi_g = dms2degrees(phi_g)
        lamb_a = dms2degrees(lamb_a)
        phi_a = dms2degrees(phi_a)
    # A multiplicação por 3600 serve para converter para segundos
    comp_meridiana = (phi_a - phi_g) * 3600
    comp_primeiro_vertical = (
        (lamb_a - lamb_g) * np.cos(np.deg2rad(phi_g))) * 3600
    return comp_meridiana, comp_primeiro_vertical


class reducoes_simplificadas(object):
    def __init__(self, DH, HA, HB, R):
        """
        @parâmetros: DH - distância horizontal
                     HA, HB - altitudes ortométricas
                     R - raio médio terrestre

        @calcula: DGC - dist. no geoide em corda
                  DGA - dist. no geoide em arco
                  DEC - dist. no elipsoide em corda
                  DEA - dist. no elipsoide em arco
        """
        self.DH = DH
        self.HA = HA
        self.HB = HB
        self.R = R

        self.DGC = (self.DH/self.R)*(self.R-(np.mean(np.abs(self.HA-self.HB))))
        self.DGA = self.DGC * (1+(self.DGC**2/(24*self.R**2)))
        self.DEC = self.DGC + (1.027*self.DGC**3*1.0E-15)
        self.DEA = self.DEC*(1+(self.DEC**2/(24*self.R**2)))


class reducoes(object):
    def __init__(self, DI, hA, hB, HA, HB, R):
        """
        @parâmetros: DI - distância inclinada
                     hA, hB - altitudes geométricas
                     HA, HB - altitudes ortométricas
                     R - raio médio terrestre

        @calcula: DGC - dist. no geoide em corda
                  DGA - dist. no geoide em arco
                  DEC - dist. no elipsoide em corda
                  DEA - dist. no elipsoide em arco
        """
        self.DI = DI
        self.hA = hA
        self.hB = hB
        self.HA = HA
        self.HB = HB
        self.R = R

        self.DEC = np.sqrt((self.DI**2-(self.hA-self.hB)**2) /
                           ((1+(self.hA/self.R))*(1+(self.hB/self.R))))
        self.DEA = self.DEC * (1+(self.DEC**2/(24*self.R**2)))
        self.DGC = np.sqrt((self.DI**2-(self.HA-self.HB)**2) /
                           ((1+(self.HA/self.R))*(1+(self.HB/self.R))))
        self.DGA = self.DGC * (1+(self.DGC**2/(24*self.R**2)))
