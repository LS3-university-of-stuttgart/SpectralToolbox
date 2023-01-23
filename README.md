# SpectralToolbox

==========================
FFT-based kriging package
==========================
Version 1.5 (01 august 2007 / WN)

--------------------------
COPYRIGHTS NOTICE
--------------------------
This software may be used for non-commercial purposes and for
university-related research and education under the condition
that the publication mentioned below will be appropriately
acknowledged and cited in reports and other publications.

Neither original nor modified versions of this code as a whole or
in parts may be passed on to others without permission of the author.

A distributable software package (e.g. under the GPL) may become
available later on.

Copyright 2007 by Wolfgang Nowak and Jochen Fritz.

--------------------------
DISCLAIMER
--------------------------
In the current version of the code, quality and consistency of input
parameters is not checked upon startup and may result in run-time
errors or erroneous results. Please read the help lines of the main
file carefully and test small problems before gaining confidence in
how to handle this code.

The authors give no warranty for the correctness of the results.

--------------------------
CONTACT INFORMATION
--------------------------
Affiliation:
  Institute of Hydraulic Engineering (IWS)
  Chair for Hydraulics and Modelling of Hydrosystems (LH2)
  University of Stuttgart, Germany
Email:
  wolfgang.nowak@iws.uni-stuttgart.de

Please report bugs to this email address.

General support and helpdesk services can not be offered.

--------------------------
INFORMATION ON METHODS
--------------------------
Detailed information on the method used here:

J. Fritz, W. Nowak and I. Neuweiler: "FFT-based Algorithms for Kriging",
submitted to Mathematical Geology (2007)

Other useful publications:

W. Nowak, S. Tenkleve and O.A. Cirpka: "Efficient computation of linearized
cross-covariance and auto-covariance matrices of interdependent quantities",
Math. Geol., 35(1), 53-66 (2003).

W. Nowak.: "Geostatistical Methods for the Identification of Flow and Transport
Parameters in Subsurface Flow", http://elib.uni-stuttgart.de/opus/frontdoor.php?source_opus=2275
(2005).

R. H. Chan and M. K. Ng: "Conjugate Gradient Methods for Toeplitz Systems",
SIAM Review, 38(3), 427-482 (1996).

J. R. Shewchuk: "An Introduction to the Conjugate Gradient Method Without
the Agonizing Pain", http://www.cs.berkeley.edu/~jrs/ (1994).

P. K. Kitanidis: "Analytical expressions of conditional mean, covariance, and
sample functions in geostatistics", Stoch. Hydrol. Hydraul., 12, 279-294 (1996).

--------------------------
CONTENTS OF THIS TOOLBOX
--------------------------
general_kriging: 	the kriging package containing all other required
			components (not the fft-based solvers) as subfunctions
general_kriging_test:	a small test routine that calls general_kriging.
			also serves as an example for how to define the input parameters
sts_pcg:		fft-based pcg solver for Toeplitz systems (regular grids)
gsts_pcg:		fft-based solvers for almost-Toeplitz systems (irregular grids)