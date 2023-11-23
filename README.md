
SpectralToolbox
==========================

A FFT-based kriging package

Version 1.5 (01 august 2007 / WN)

--------------------------
COPYRIGHTS NOTICE
--------------------------

See [License](https://github.com/LS3-university-of-stuttgart/SpectralToolbox/blob/main/LICENSE.md)

Please cite the following paper when using this software.

```biblatex
@article{fritz2009,
	title        = {Application of FFT-based Algorithms for Large-Scale Universal Kriging Problems},
	author       = {Fritz, J. and Neuweiler, I. and Nowak, W.},
	year         = 2009,
	month        = {Jul},
	day          = {01},
	journal      = {Mathematical Geosciences},
	volume       = 41,
	number       = 5,
	pages        = {509--533},
	doi          = {10.1007/s11004-009-9220-x},
	issn         = {1874-8953},
	url          = {https://doi.org/10.1007/s11004-009-9220-x},
	abstract     = {Looking at kriging problems with huge numbers of estimation points and measurements, computational power and storage capacities often pose heavy limitations to the maximum manageable problem size. In the past, a list of FFT-based algorithms for matrix operations have been developed. They allow extremely fast convolution, superposition and inversion of covariance matrices under certain conditions. If adequately used in kriging problems, these algorithms lead to drastic speedup and reductions in storage requirements without changing the kriging estimator. However, they require second-order stationary covariance functions, estimation on regular grids, and the measurements must also form a regular grid. In this study, we show how to alleviate these rather heavy and many times unrealistic restrictions. Stationarity can be generalized to intrinsicity and beyond, if decomposing kriging problems into the sum of a stationary problem and a formally decoupled regression task. We use universal kriging, because it covers arbitrary forms of unknown drift and all cases of generalized covariance functions. Even more general, we use an extension to uncertain rather than unknown drift coefficients. The sampling locations may now be irregular, but must form a subset of the estimation grid. Finally, we present asymptotically exact but fast approximations to the estimation variance and point out application to conditional simulation, cokriging and sequential kriging. The drastic gain in computational and storage efficiency is demonstrated in test cases. Especially high-resolution and data-rich fields such as rainfall interpolation from radar measurements or seismic or other geophysical inversion can benefit from these improvements.}
}
```

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
  Stochastic Simulation and Safety Research for Hydrosystems (LS3)
  University of Stuttgart, Germany

[Prof. Nowak](https://www.iws.uni-stuttgart.de/institut/team/Nowak-00003/)

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
