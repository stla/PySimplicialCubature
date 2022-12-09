from pysimplicialcubature.simplicialcubature import integratePolynomialOnSimplex, integrateOnSimplex
from sympy import Poly
from sympy.abc import x, y, z

def test_polynomial():
    # simplex vertices
    v1 = [1.0, 1.0, 1.0] 
    v2 = [2.0, 2.0, 3.0] 
    v3 = [3.0, 4.0, 5.0] 
    v4 = [3.0, 2.0, 1.0]
    # simplex
    S = [v1, v2, v3, v4]
    # polynomial to integrate
    P = Poly(x**4 + y + 2*x*y**2 - 3*z, x, y, z, domain = "RR")
    # integral of P
    I_P = integratePolynomialOnSimplex(P, S)
    # polynomial P as function
    f = lambda x : x[0]**4 + x[1] + 2*x[0]*x[1]**2 - 3*x[2]
    # integral of f
    I_f = integrateOnSimplex(f, S)
    # compare results
    assert abs(I_P - I_f["integral"]) < 1e-5
