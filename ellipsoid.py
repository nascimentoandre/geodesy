from numpy import sin, sqrt, radians


class Ellipsoid(object):
    """
    A class used to represent an ellipsoid.
    For initialization, it requires:
    a -> major semi-axis
    f -> flattening
    """

    def __init__(self, a, f):
        """
        Parameters
        ----------
        a : float
            Major semi-axis
        f : float
            Flattening
        b : float
            Minor semi-axis
        e1 : float
            First eccentricity
        e2 : float
            Second eccentricity
        ---------
        """
        self.a = a
        self.f = f
        self.b = self.a - self.f*self.a
        self.e1 = sqrt((self.a**2-self.b**2)/self.a**2)
        self.e2 = sqrt((self.a**2-self.b**2)/self.b**2)

    def medirianNormal(self, phi):
        """
        Calculates the normal to point P with latitude phi on ellipsoid surface,
        from the meridian plane

        Parameters
        ----------
        phi : float
            Latitude in degrees
        """
        mN = self.a/(sqrt(1-(self.e1**2)*sin(radians(phi))**2))
        return mN

    def equatorNormal(self, phi):
        """
        Calculates the normal to point P with latitude phi on ellipsoid surface,
        from the equator plane

        Parameters
        ----------
        phi : float
            Latitude in degrees
        """
        eN = self.medirianNormal(phi)*(1-self.e1**2)
        return eN


def export_ellipsoids():
    """
    Returns a dictionary of ellipsoids from the ellipsoids.txt file
    """

    with open("ellipsoid.txt") as file:
        ellipsoid_dict = {}
        for line in file:
            a = float(line.split()[0])
            f_str_inv = line.split()[1].split("/")[1]
            f = 1/float(f_str_inv)
            ellipsoid_dict[line.split()[-1]] = Ellipsoid(a, f)
    return ellipsoid_dict
