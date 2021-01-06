import numpy as np


def astrloc2geodloc(X0, Y0, Z0, Xp, Yp, Zp, lamb0, phi0):
    """
    astrloc2geodloc(X0, Y0, Z0, Xp, Yp, Zp, lamb0, phi0)

    Permite passar coordenadas do sistema astronômico local (sistema topográfico) para o sistema geodésico
    local (sistema topocêntrico).

    Parâmetros
    --------------
    X0: Coordenada topográfica X do ponto datum.
    Y0: Coordenada topográfica Y do ponto datum.
    Z0: Coordenada topográfica Z do ponto datum.
    Xp: Coordenada topográfica X do ponto considerado.
    Yp: Coordenada topográfica Y do ponto considerado.
    Zp: Coordenada topográfica Z do ponto considerado.
    lamb0: Longitude geodésica do ponto datum.
    phi0: Latitude geodésica do ponto datum.

    Retorna
    ---------------
    ep: Coordenada topocêntrica leste do ponto considerado.
    np_: Coordenada topocêntrica norte do ponto considerado.
    up: Coordenada topocêntrica vertical do ponto considerado.
    """
    lamb0 = np.deg2rad(lamb0)
    phi0 = np.deg2rad(phi0)
    v_astrloc = np.array([[Xp-X0], [Yp-Y0], [Zp-Z0]])
    matriz_transf = np.array([[-np.sin(lamb0), np.cos(lamb0), 0],
                              [-np.sin(phi0)*np.cos(lamb0), -np.sin(phi0)
                               * np.sin(lamb0), np.cos(phi0)],
                              [np.cos(phi0)*np.cos(lamb0), np.cos(phi0)*np.sin(lamb0), np.sin(phi0)]])
    ep, np_, up = np.matmul(matriz_transf, v_astrloc)
    return ep[0], np_[0], up[0]


def geodloc2astrloc(e, n, u, lamb0, phi0):
    """
    geodloc2astrloc(e, n, u, lamb0, phi0)

    Permite passar coordenadas do sistema geodésico local (sistema topocêntrico) para o sistema astronômico
    local (sistema topográfico).

    Parâmetros
    --------------
    e: Variação da coordenada topocêntrica leste.
    n: Variação da coordenada topocêntrica norte.
    u: Variação da coordenada topocêntrica vertical.
    lamb0: Longitude geodésica do ponto datum.
    phi0: Latitude geodésica do ponto datum.

    Retorna
    ---------------
    dx: Variação da coordenada topográfica X.
    dy: Variação da coordenada topográfica Y.
    dz: Variação da coordenada topográfica Z.
    """
    lamb0 = np.deg2rad(lamb0)
    phi0 = np.deg2rad(phi0)
    v_geodloc = np.array([[e], [n], [u]])
    matriz_transf = np.array([[-np.sin(lamb0), -np.sin(phi0)*np.cos(lamb0), np.cos(phi0)*np.cos(lamb0)],
                              [np.cos(lamb0), -np.sin(phi0)*np.sin(lamb0),
                               np.cos(phi0)*np.sin(lamb0)],
                              [0, np.cos(phi0), np.sin(phi0)]])
    dx, dy, dz = np.matmul(matriz_transf, v_geodloc)
    return dx[0], dy[0], dz[0]
