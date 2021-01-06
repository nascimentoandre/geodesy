from coordinate_conv import geod2cart, cart2geod

# Criando um dicionário com dicionários de parâmetros de transformação entre sistemas de referência

parametros = {'SIRGAS20002SAD69': {'dx': 67.35, 'dy': -3.88, 'dz': 38.22},
              'SAD692SIRGAS2000': {'dx': -67.35, 'dy': +3.88, 'dz': -38.22},
              'SICAD2SIRGAS2000': {'dx': -144.35, 'dy': 242.88, 'dz': -33.22},
              'SIRGAS20002SICAD': {'dx': 144.35, 'dy': -242.88, 'dz': 33.22}}


def conv_sist_ref(lamb, phi, h, elip1, elip2, dms=False):
    """
    conv_sist_ref(lamb, phi, h, elip1, elip2)

    Permite converter coordenadas entre os diferentes sistemas de referência definidos. Sistemas de
    referência definidos até o momento:
    1) SIRGAS2000
    2) SAD69
    3) SICAD

    Parâmetros
    -----------
    lamb: Valor da longitude a ser convertida.
    phi: Valor da latitude a ser convertida.
    h: Valor da altitude geométrica a ser convertida.
    elip1: Sistema de referência no qual as coordenadas estão.
    elip2: Sistema de referência para o qual se deseja converter as coordenadas.
    dms: Formato das coordenadas de entrada. Caso as coordenadas estejam em graus, minutos e segundos,
    deve-se atribuir um valor True.

    Retorna
    ---------
    lamb: Valor da longitude no novo sistema de referência.
    phi: Valor da latitude no novo sistema de referência.
    h: Valor da altitude geométrica no novo sistema de referência.
    """
    global parametros
    if dms:
        X, Y, Z = geod2cart(lamb, phi, h, elip1, dms=True)
    else:
        X, Y, Z = geod2cart(lamb, phi, h, elip1)
    modo = elip1.nome + '2' + elip2.nome
    X2, Y2, Z2 = X + parametros[modo]['dx'], Y + \
        parametros[modo]['dy'], Z + parametros[modo]['dz']
    if dms:
        return cart2geod(X2, Y2, Z2, elip2, dms=True)
    else:
        return cart2geod(X2, Y2, Z2, elip2)
