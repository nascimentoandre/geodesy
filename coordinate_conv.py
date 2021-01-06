import numpy as np


def dms2degrees(coord):
    """
    dms2degrees(coord)

    Converte uma coordenada em graus, minutos e segundos para graus decimais.

    Parâmetros
    -----------
    coord: Lista (ou tupla) com os valores dos graus, minutos e segundos. Exemplo: [42, 55, 1.9845].

    Retorna
    ---------
    degrees: Valor da coordenada em graus decimais. 
    """
    if coord[0] < 0:
        degrees = coord[0] - coord[1]/60 - coord[2]/3600
    else:
        degrees = coord[0] + coord[1]/60 + coord[2]/3600
    return degrees


def degrees2dms(coord):
    """
    degrees2dms(coord)

    Converte uma coordenada em graus decimais para graus, minutos e segundos.

    Parâmetros
    -----------
    coord: Valor da coordenada em graus decimais.

    Retorna
    ---------
    Lista com os valores dos graus, minutos e segundos.
    """
    coordstr = str(coord).split('.')[1]
    resto1 = float('0.'+coordstr)
    m = resto1 * 60
    minutos = float(str(m).split('.')[0])
    resto2 = float('0.'+str(m).split('.')[1])
    segundos = resto2 * 60
    return [float(str(coord).split('.')[0]), minutos, segundos]


def geod2cart(lamb, phi, h, elip, dms=False):
    """
    geod2cart(lamb, phi, h, elip, dms=False)

    Converte coordenadas geodésicas em coordenadas cartesianas tridimensionais.

    Parâmetros
    -----------
    lamb: Valor da longitude.
    phi: Valor da latitude.
    h: Valor da altitude geométrica.
    elip: Qual o sistema de refência em questão.
    dms: Atributo relacionado ao formato da latitude e da longitude. O padrão é False, o que significa
    que as coordenadas estão em graus decimais. Caso elas estejam em graus, minutos e segundos, é 
    necessário atribuir o valor True.

    Retorna
    ---------
    X: Valor da coordenada cartesiana na direção X.
    Y: Valor da coordenada cartesiana na direção Y.
    Z: Valor da coordenada cartesiana na direção Z.
    """
    # Caso a coordenada esteja em graus, minutos e segundos, precisaremos convertê-la para graus decimais
    if dms:
        lamb = dms2degrees(lamb)
        phi = dms2degrees(phi)
    X = (elip.grandeNormal(phi) + h) * \
        np.cos(np.radians(phi)) * np.cos(np.radians(lamb))
    Y = (elip.grandeNormal(phi) + h) * \
        np.cos(np.radians(phi)) * np.sin(np.radians(lamb))
    Z = (elip.grandeNormal(phi)*(1 - elip.e1**2) + h) * np.sin(np.radians(phi))
    return X, Y, Z


def cart2geod(X, Y, Z, elip, dms=False):
    """
    cart2geod(X, Y, Z, elip, dms=False)

    Converte coordenadas cartesianas tridimensionais em coordenadas geodésicas.

    Parâmetros
    -----------
    X: Valor da coordenada cartesiana na direção X.
    Y: Valor da coordenada cartesiana na direção Y.
    Z: Valor da coordenada cartesiana na direção Z.
    elip: Qual o sistema de refência em questão.
    dms: Atributo relacionado ao formato da saída. O padrão é False, o que significa que as
    coordenadas irão sair em graus decimais. Caso o usuário deseje a saída em graus, minutos
    e segundos, é necessário atribuir o valor True.

    Retorna
    ---------
    lamb: Valor da longitude.
    phi: Valor da latitude.
    h: Valor da altitude geométrica.
    """
    u = np.arctan((Z/(np.sqrt(X**2+Y**2)))*(elip.a/elip.b))
    lamb = np.rad2deg(np.arctan(Y/X))
    phi = np.rad2deg(np.arctan((Z + elip.e2**2*elip.b*np.sin(u)**3)/((np.sqrt(X**2+Y**2)
                                                                      - elip.e1**2*elip.a*np.cos(u)**3))))
    h = (np.sqrt(X**2+Y**2)/np.cos(np.deg2rad(phi))) - elip.grandeNormal(phi)
    # A verificação abaixo faz a análise de quadrantes e corrige o valor de lambda
    if X < 0 and Y < 0:
        lamb = -(180-lamb)
    elif X < 0 and Y >= 0:
        lamb = 180 + lamb
    # Caso deseje que a saída esteja em graus, minutos e segundos, basta atribuir o valor True ao atributo dms
    if dms:
        lamb = degrees2dms(lamb)
        phi = degrees2dms(phi)
    return lamb, phi, h
