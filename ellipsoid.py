from math import sin, sqrt, radians


class Elipsoide(object):
    """
    Esta classe cria um objeto da classe elipsoide
    """

    def __init__(self, a, f):
        """
        O método construtor recebe os valores do semieixo maior (a) e do achatamento (f) do elipsoide,
        e a partir deles calcula os valores do semieixo menor (b), da primeira excentricidade (e1) e
        da segunda excentricidade (e2). Além disso, um dos atributos do sistema de referência a ser informado
        é o seu nome, que será usado unicamente para facilitar a conversão de coordenadas entre diferentes
        sistemas de referência.
        """
        self.a = a
        self.f = f
        self.b = self.a - self.f*self.a
        self.e1 = sqrt((self.a**2-self.b**2)/self.a**2)
        self.e2 = sqrt((self.a**2-self.b**2)/self.b**2)

    def grandeNormal(self, phi):
        """
        Calcula a grande normal ao elipsoide para um determinado valor de latitude phi (em graus)
        """
        gN = self.a/(sqrt(1-(self.e1**2)*sin(radians(phi))**2))
        return gN

    def pequenaNormal(self, phi):
        """
        Calcula a pequena normal ao elipsoide para um determinado valor de latitude phi (em graus)
        """
        pN = self.grandeNormal(phi)*(1-self.e1**2)
        return pN


def export_ellipsoids():
    with open("ellipsoid.txt") as file:
        ellipsoid_dict = {}
        for line in file:
            a = float(line.split()[0])
            f_str_inv = line.split()[1].split("/")[1]
            f = 1/float(f_str_inv)
            ellipsoid_dict[line.split()[-1]] = Elipsoide(a, f)
    return ellipsoid_dict
