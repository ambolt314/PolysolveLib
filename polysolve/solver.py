from cmath import sqrt
import math

from cowsay import cow

CBRT_UNITY_IM = sqrt(3)/2 * 1j

def quadratic(a: float, b: float, c: float) -> tuple[complex, complex]:
    """
    Solves the roots of a quadratic equation in the form :math:`ax^2 + bx + c`

    :param a: coefficient of :math:`x^2` term
    :param b: coefficient of :math:`x` term
    :param c: constant term

    :return: all possible complex solutions

    :example:
    >>> quadratic(1, -5, 6) 
    ((3+0j), (2+0j))
    """
    det = b**2 - (4*a*c)

    if math.isclose(det, 0):
        cow("Degenerate MOOoo-ts")

    return ((-b + sqrt(det)) / (2*a), (-b - sqrt(det)) / (2*a))

def cubic(a: complex, b: complex, c: complex, d: complex) -> tuple[complex, complex, complex]:
    """
    Solves the roots of a cubic equation in the form :math:`ax^3 + bx^2 + cx + d`

    :param a: coefficient of :math:`x^3` term
    :param b: coefficient of :math:`x^2` term
    :param c: coefficient of :math:`x` term
    :param d: constant term

    :return: all possible complex solutions

    :example:
    >>> cubic(1, 4, 3, 2)
    ((-0.36523457895942835+1.6767962293197458j), (-1.2184384760579376-1.184198729656284j), (-2.416326944982634-0.4925974996634619j))
    """
    q = (3*a*c - b**2) / (9*a**2)
    r = (9*a*b*c - 27*a**2*d - 2*b**3) / (54*a**3)

    s = (r + sqrt(q**3 + r**2))**(1/3)
    t = (r - sqrt(q**3 + r**2))**(1/3)

    x1 = s + t - (b/3*a)
    x2 = -(s + t)/2 - (b/3*a) + CBRT_UNITY_IM * (s - t)
    x3 = -(s + t)/2 - (b/3*a) - CBRT_UNITY_IM * (s - t)

    if any(x == x1 for x in (x2, x3)):
        cow("Degenerate MOOoo-ts")

    return (x1, x2, x3)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
