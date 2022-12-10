.. pysimplicialcubature documentation master file, created by
   sphinx-quickstart on Fri Dec  9 20:41:25 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pysimplicialcubature's documentation!
================================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

This package is a port of a part of the R package **SimplicalCubature**, 
written by John P. Nolan, and which contains R translations of 
some Matlab and Fortran code written by Alan Genz. In addition it 
provides a function for the exact computation of the integral of a 
polynomial over a simplex.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


Members functions
=================

.. automodule:: pysimplicialcubature.simplicialcubature
    :members:


References
==========

* A. Genz and R. Cools. 
*An adaptive numerical cubature algorithm for simplices.* 
ACM Trans. Math. Software 29, 297-308 (2003).

* Jean B. Lasserre.
*Simple formula for the integration of polynomials on a simplex.* 
BIT Numerical Mathematics 61, 523-533 (2021).