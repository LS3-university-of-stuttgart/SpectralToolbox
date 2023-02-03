clear all
clc
clf

% setus up input data for GENERAL_KRIGING and runs a small test problem
% (should take a few seconds on a machine purchased in 2007)

% definition of prior covariance
s.model      = 'gaussian';                         % geostatistical model for unknowns
s.variance   = 1;                                  % geostatistical parameter for exponential model
s.lambda     = [6 6 2];                            % geostatistical parameter for exponential model
s.nugget     = 0;                                  % nugget effect
s.micro      = 0;                                  % microscale smoothing parameter (before nugget)

% definition of the grid for the unknowns s
s.n_pts      = [100 100 50 ];                      % number of unknowns in each direction
s.d_pts      = [0.4 0.4 0.2];                      % grid spacing in each direction
s.npts       = prod(s.n_pts);                      % number of unknowns (not required unless used in problem generation)

% definition of mean
b.model      = 'known';                            % uncertain, known, zero or unknown
b.n          = 1;                                  % number of base functions
b.beta_pri   = 0        ;                          % prior mean coefficients for trend functions
b.Qbb        = 0;                                  % uncertainty of prior mean: covariance matrix for trend coefficients
b.function   = ones(s.npts,1);                     % base functions (vector of ones for constant mean)

% definition of measurement locations for measurements y
y.gridtype   = 'irregular';                        % type of measurement grid (regular or irregular)
y.npts       = 50;                                 % number of observations (not required unless used for problem generation)
[zzz,aux]    = sort(randn(s.npts,1));              % randomized choice of locations
y.indices    = aux(1:y.npts);                      % measurement indices in field of unknowns (required for irregular grids)
% y.indices    = round(s.npts/2)-60-120*50+(120*100)*[-24:25]';  % measurement indices in field of unknowns (required for irregular grids)
% y.gridtype   = 'regular';                          % type of measurement grid (regular or irregular)
% y.d_ratio    = [10 10 10];                         % resampling ratio (required for regular grids)
% y.x_first    = [5  5  5 ];                         % location of first measurement: grid node number of s (required for regular grids only)
% y.n_pts      = [5  5  5 ];                         % size of subsampled field (required for regular grids)
% y.npts       = prod(y.n_pts);                      % number of observations (not required unless used for problem generation)

% generation of artificial data set y
% y.error      = 1.^2;                               % measurement error (scalar) expressed as variance
% y.values     = randn(y.npts,1)*sqrt(s.variance)...
%                + b.beta_pri;                       % randomized measurement values
y.error      = 0.0001.^2;                               % measurement error (scalar) expressed as variance
y.values     = ones(y.npts,1) + b.beta_pri;           % randomized measurement values

% kriging method options
options.superpos = 'fft';                          % superposition method: fft or standard
options.solver   = 'standard';                     % solver method: fft (fft-reg, fft-irreg) or standard
options.estvar   = 'none';                         % estimation variance method: full, one-point, speedy or none
options.plot     = true;                           % plotting flag: true or false

% specific options for fft-based solvers
options.tol      = 1e-10;                          % solver relative residual (required for fft-solvers: usually about 1e-10)
options.maxit    = 200;                            % solver maximum iteration number (required for fft-solvers: usually number of measurements or less)
options.cond     = 1e6;                            % solver regularization parameter (required for fft-solvers: usually 1e6)
options.verbose  = 0;                              % solver verbosity (required for fft-solvers: 0, 1 or 2)
options.kalstr   = 0;                              % solver filter strength (required for fft-solvers: zero or some percent of variance)
options.flag_Strang = 0;                           % solver method flag (required for fft-solvers: 1 or 0)
options.maxprime = 7;                              % embedding optimization parameter (2,3,5,7,...)

[estimate,ksi,beta,est_var]=general_kriging(s,b,y,options);