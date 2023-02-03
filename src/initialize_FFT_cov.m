function [FFTQe,se] = initialize_FFT_cov(s,options)
% INITIALIZE_FFT_COV initializes the variable FFTQe for later recycling
% FFTQe is the FFTN of the embedded version of the first row
% of a stationary covariance matrix (symmetric block Toeplitz structure)
% version 29 january 2008 / WN
%
% Required geostatistical input parameters:
% s.name           : parametric geostatistical model name (string code)
% s.variance       : variance of model (nonnegative scalar)
% s.lambda         : vector of correlation length in y,x,z-direction
%
% Required geometrical input parameters:
% s.d_pts          : length of grid points in y,x,z-direction
% s.n_pts          : number of grid points in y,x,z-direction
%
% Optional input paramters:
% s.micro          : microscale smoothing parameter relative to lambda
%                    (scalar, after Kitanidis book). default value is
%                    zero
% s.nugget         : nugget effect (additional, variance) on top of chosen model (scalar)
% s.kappa          : shape parameter for Matern covariance model.
%                    default value is 1.5
% s.periodicity    : vector of flags setting periodicity in
%                    y,x,z-direction. default value is [0 0 0]
%
% output parameters:
% FFTQe


% checking input parameter
s = check_structure_s(s);

% finding appropriate size for the embedding (approximately
% ensuring positive definiteness of embedded covariance matrix)
se     = find_embedding(s,options);

% get the effective (anisotropic) separation distances
h_effe = evaluate_separation(se.x_pts,s.lambda,s.rotation,s.micro,s.nd);

% get the covariance function (embedded version)
Qe     = evaluate_covariance(s,h_effe);

% computing FFTn and eigenvalues of embedded covariance matrix
FFTQe          = abs(fftn(reshape(Qe,[se.n_pts,1])));