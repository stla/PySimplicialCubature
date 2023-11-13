# pysimplicialcubature

<!-- badges: start -->
[![Documentation status](https://readthedocs.org/projects/pysimplicialcubature/badge/)](http://pysimplicialcubature.readthedocs.io)
<!-- badges: end -->

This package is a port of a part of the R package **SimplicialCubature**, 
written by John P. Nolan, and which contains R translations of 
some Matlab and Fortran code written by Alan Genz. In addition it 
provides a function for the exact computation of the integral of a 
polynomial over a simplex.

___

A simplex is a triangle in dimension 2, a tetrahedron in dimension 3. 
This package provides two main functions: `integrateOnSimplex`, to integrate 
an arbitrary function on a simplex, and `integratePolynomialOnSimplex`, to 
get the exact value of the integral of a multivariate polynomial on a 
simplex.

Suppose for example you want to evaluate the following integral:

$$\int\_0^1\int\_0^x\int\_0^y \exp(x + y + z) \text{d}z \text{d}y \text{d}x.$$

```python
from pysimplicialcubature.simplicialcubature import integrateOnSimplex
from math import exp

# simplex vertices
v1 = [0.0, 0.0, 0.0] 
v2 = [1.0, 1.0, 1.0] 
v3 = [0.0, 1.0, 1.0] 
v4 = [0.0, 0.0, 1.0]
# simplex
S = [v1, v2, v3, v4]
# function to integrate
f = lambda x : exp(x[0] + x[1] + x[2])
# integral of f on S
I_f = integrateOnSimplex(f, S)
I_f["integral"]
# 0.8455356728324119
```

The exact value of this integral is ${(e-1)}^3/6 \approx 0.8455356852954753$.

Now let's turn to a polynomial example. You have to define the polynomial with 
`sympy.Poly`.

```python
from pysimplicialcubature.simplicialcubature import integratePolynomialOnSimplex
from sympy import Poly
from sympy.abc import x, y, z

# simplex vertices
v1 = [1.0, 1.0, 1.0] 
v2 = [2.0, 2.0, 3.0] 
v3 = [3.0, 4.0, 5.0] 
v4 = [3.0, 2.0, 1.0]
# simplex
S = [v1, v2, v3, v4]
# polynomial to integrate
P = Poly(x**4 + y + 2*x*y**2 - 3*z, x, y, z, domain = "RR")
# integral of P on S
integratePolynomialOnSimplex(P, S)
# -0.253571428571429
```

We can do better: by defining the vertex coordinates as fractions, and by 
taking the field of rational numbers as the domain of the polynomial, we will 
get the exact value of the integral, without numerical error.

```python
from pysimplicialcubature.simplicialcubature import integratePolynomialOnSimplex
from sympy import Poly
from sympy.abc import x, y, z
from fractions import Fraction

# simplex vertices
v1 = [Fraction(0), Fraction(0), Fraction(0)] 
v2 = [Fraction(1), Fraction(1), Fraction(1)] 
v3 = [Fraction(0), Fraction(1), Fraction(1)] 
v4 = [Fraction(0), Fraction(0), Fraction(1)]
# simplex
S = [v1, v2, v3, v4]
# polynomial to integrate
P = Poly(x**4 + y + 2*x*y**2 - 3*z, x, y, z, domain = "RR")
# integral of P on S
integratePolynomialOnSimplex(P, S)
# -71/280 (= -0.25357142857142856...)
```


## References

- A. Genz and R. Cools. 
*An adaptive numerical cubature algorithm for simplices.* 
ACM Trans. Math. Software 29, 297-308 (2003).

- Jean B. Lasserre.
*Simple formula for the integration of polynomials on a simplex.* 
BIT Numerical Mathematics 61, 523-533 (2021).