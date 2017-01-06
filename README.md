Locally one dimensionnal heat transfer.
=======================================

This repository contains the Matlab code to solve a 2D heat transfer problem
using the locally one dimensional (or alternate direction implicit) method
described in `Levy, Hoang & Le Corre, Materials Sciences and Applications, 
2017, 8`.
This handles temperature dependent material properties (ie nonlinear problem)
but resolution is done in linear since time steps are small such that
non-linearity is limited
$P_z$ is a series a 1D problem and $P_m$ one single 1D problem.

The algorithm is described in figure algorithm.gif.