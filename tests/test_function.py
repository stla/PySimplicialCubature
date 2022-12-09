from pysimplicialcubature.simplicialcubature import integrateOnSimplex
from math import exp

def test_function():
    # simplex vertices
    v1 = [0.0, 0.0, 0.0] 
    v2 = [1.0, 1.0, 1.0] 
    v3 = [0.0, 1.0, 1.0] 
    v4 = [0.0, 0.0, 1.0]
    # simplex
    S = [v1, v2, v3, v4]
    # function to integrate
    f = lambda x : exp(x[0] + x[1] + x[2])
    # integral of f
    I_f = integrateOnSimplex(f, S)
    # compare results
    assert abs((exp(1) - 1)**3 / 6 - I_f["integral"]) < 1e-5
