Locally one dimensionnal heat transfer.
=======================================

This repository contains the Matlab code to solve a 2D heat transfer problem
using the locally one dimensional (or alternate direction implicit) method
described in `Hoang, Levy & Le Corre, International Journal for Material Forming, 
2016`.
This handles temperature dependent material properties (ie nonlinear problem)
$P_z$ is a series a 1D problem and $P_m$ one single 1D problem.

The algorithm is described in figure algorithm.gif.