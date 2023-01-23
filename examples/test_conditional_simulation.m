addpath('../src')
% setus up input data for GENERAL_KRIGING and runs a small test problem
% (should take a few seconds on a machine purchased in 2007)

% definition of prior covariance
s.model      = 'gaussian';                         % geostatistical model for unknowns
s.variance   = 4;                                  % geostatistical parameter for exponential model
s.lambda     = [10  10  1];                        % geostatistical parameter for exponential model
s.nugget     = 0;                                  % nugget effect
s.kappa      = 1.5;                                % shape parameter (for matern covariance only)
s.micro      = 0;                                  % microscale smoothing parameter (before nugget)

% definition of the grid for the unknowns s
s.n_pts      = [350 350 4];                       % number of unknowns in each direction
s.d_pts      = [1   1   1];                        % grid spacing in each direction
s.npts       = prod(s.n_pts);                      % number of unknowns (not required unless used in problem generation)
s            = ndgrid_setup(s);                    %%just so that you have the full grid info available

% definition of mean
b.model      = 'known';                            % uncertain, known, zero or unknown
b.n          = 1;                                  % number of base functions
b.beta_pri   = 0;                                  % prior mean coefficients for trend functions
b.Qbb        = 0;                                  % uncertainty of prior mean: covariance matrix for trend coefficients
b.trend{1}   = ones(s.n_pts);                      % constant mean function
b.trends     = b.trend{1}(:);

% definition of measurement locations for measurements y
y.gridtype   = 'irregular';                        % type of measurement grid (regular or irregular)
y.npts       = 10;                                 % number of observations (not required unless used for problem generation)
[zzz,aux]    = sort(randn(s.npts,1));              % randomized choice of locations
y.indices    = aux(1:y.npts);                      % measurement indices in field of unknowns (required for irregular grids)
% y.gridtype   = 'regular';                          % type of measurement grid (regular or irregular)
% y.d_ratio    = [10 10];                            % resampling ratio (required for regular grids)
% y.x_first    = [10 10];                            % location of first measurement: grid node number of s (required for regular grids only)
% y.n_pts      = [10 10];                            % size of subsampled field (required for regular grids)
% y.npts       = prod(y.n_pts);                      % number of observations (not required unless used for problem generation)

% generation of artificial data set y
y.error      = 0.001.^2;                               % measurement error (scalar) expressed as variance
y.values     = 3*ones(y.npts,1);

% kriging method options
options.superpos = 'fft';                          % superposition method: fft or standard
options.solver   = 'standard';                     % solver method: fft (fft-reg, fft-irreg) or standard
options.estvar   = 'none';                         % estimation variance method: full, one-point, speedy or none
options.plot     = false;                          % plotting flag: true or false

% specific options for fft-based solvers
options.tol      = 1e-10;                          % solver relative residual (required for fft-solvers: usually about 1e-10)
options.maxit    = 200;                            % solver maximum iteration number (required for fft-solvers: usually number of measurements or less)
options.cond     = 1e6;                            % solver regularization parameter (required for fft-solvers: usually 1e6)
options.verbose  = 0;                              % solver verbosity (required for fft-solvers: 0, 1 or 2)
options.kalstr   = 0;                              % solver filter strength (required for fft-solvers: zero or some percent of variance)
options.flag_Strang = 0;                           % solver method flag (required for fft-solvers: 1 or 0)
options.maxprime = 7;                              % embedding optimization parameter (2,3,5,7,...)

% generate a random realization of the field as true field
Y_true           = generate_randomfield(s,b,options);

% generate a random realization of the field
Y                = generate_randomfield(s,b,options);
s                = ndgrid_setup(s); %%just so that you have the full grid info available

% generate the (randomized) residuals at the kriging points
y.values         = y.values-Y(y.indices) + randn(y.npts,1)*sqrt(y.error);

% set kriging mode to zero mean
b.beta_pri   = 0;                                  % prior mean coefficients for trend functions

[estimate,ksi,beta,est_var,s] = general_kriging(s,b,y,options);

plotter_nd(s.x_vec{1},s.x_vec{2},s.x_vec{3},Y+estimate,'conditional sim',[1 1 1],s.n_pts,[],[],0,100,0);
hold on
plot(s.x_pts{2}(y.indices),s.x_pts{1}(y.indices),'xk','markersize',10)
plot(s.x_pts{2}(y.indices),s.x_pts{1}(y.indices),'ok','markersize',10)
hold off
caxis([-3 3])