from pysimplicialcubature.simplicialcubature import integratePolynomialOnSimplex
from sympy import Poly
from sympy.abc import x, y, z
from fractions import Fraction

def test_polynomialQQ():
    # simplex vertices
    v1 = [Fraction(0), Fraction(0), Fraction(0)] 
    v2 = [Fraction(1), Fraction(1), Fraction(1)] 
    v3 = [Fraction(0), Fraction(1), Fraction(1)] 
    v4 = [Fraction(0), Fraction(0), Fraction(1)]
    # simplex
    S = [v1, v2, v3, v4]
    # polynomial to integrate
    P = Poly(x**4 + y + 2*x*y**2 - 3*z, x, y, z, domain = "QQ")
    # integral of P
    I_P = integratePolynomialOnSimplex(P, S)
    # assert result
    assert I_P == Fraction(-71, 280)
