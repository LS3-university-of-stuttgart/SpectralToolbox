function [estimate,ksi,beta,est_var,s] = general_kriging(s,b,y,options)
% GENERAL_KRIGING computes kriging estimate for generalized case of
% uncertain mean, using the function extimate form and a user-defined
% set of standard or fft-based methods.
%
% [ESTIMATE,KSI,BETA,EST_VAR] = GENERAL_KRIGING(S,B,Y,OPTIONS)
% returns the kriging estimate for the problem specified by the input
% parameters S, B, Y, and OPTIONS. ESTIMATE is kriging result,
% KSI and BETA are the kriging weights and EST_VAR is the
% corresponding estimation variance.
%
% Input parameters are structures for:
% S:           unknowns
% B:           base functinos
% Y:           measurements
% OPTIONS:     general options describing the problem setup and user's
%              choice of method
% 
% In the current code, quality of input parameters is not yet checked
% upon startup and may result in run-time errors later on. Please read
% the following help lines carefully and test small problems before
% gaining confidence in how to handle this code. The authors give no
% warranty for the correctness of the results.
%
% Definition of prior covariance: Subfields of structure S:
% .model       geostatistical model type:
%               gaussian, exponential, spherical (to be extended)
% .variance    geostatistical model parameter for field variance:
%               positive scalar
% .lambda      geostatistical model parameter for correlation length:
%               (may be anisotropic, rotated anisotropy not implemented)
%               positive vector of length d (d is dimensionality)
% .nugget      geostatistical model parameter for nugget effect:
%               positive scalar (adds to variance)
% .micro       microscale smoothing parameter (see book by Kitanidis):
%               positive scalar, relative to .lambda
%               (applied before adding nugget)
% 
% Definition of the grid of unknowns: Subfields of structure S:
% .n_pts       number of unknowns in each direction:
%               positive integer vector of length d
% .d_pts       grid spacing in each direction:
%               positive integer vector of length d
% 
% Definition of mean: Subfields of structure B:
% .model       model type for the stochastic mean of the unknowns:
%               uncertain, known, zero or unknown (string)
% .n           model complexity (number of base functions):
%               positive scalar
% .beta_pri    prior mean coefficients for trend functions:
%               vector of length as specified in b.n
% .Qbb         uncertainty of prior mean
%               (covariance matrix for trend coefficients):
%               positive-definite square matrix sized b.n times b.n
% .function    base functions (vector of ones for constant mean):
%               matrix sized prod(n_pts) times b.n
% 
% Definition of measurement locations: Subfields of structure Y:
% .gridtype    type of measurement grid:
%               regular or irregular (string)
% .indices     location of measurements in grid of unknowns:
%               vector with length equal to number of measurements
%               (single index notation)
%               (required for irregular grids only)
% .d_ratio     resampling ratio:
%               positive integer vector of length d
%               (required for definition of regular grids only)
% .x_first     location of first measurement in grid of unknowns:
%               positive integer vector of length d
%               (multiple subscript notation)
%               (required for regular grids only)
% .n_pts       size of subsampled field in numbers of measurements
%               (because regular grid of measurements may be smaller than domain):
%               positive integer vector of length d
%               (required for regular grids only)
% 
% data set of measurements: Subfields of structure Y:
% .error       measurement error (scalar) expressed as variance:
%               positive scalar
% .values      measurement values:
%               vector with length equal to number of measurements
% 
% kriging method options: Subfields of structure OPTIONS:
% .superpos    superposition method:
%               fft or standard (string)
% .solver      solver method:
%               fft (fft-reg, fft-irreg) or standard (string)
% .estvar      estimation variance method:
%               full, one-point, speedy or none (string) ("hybrid" may be
%               added later)
% .maxprime    embedding optimization parameter:
%               2,3,5,7,... (required for both standard and fft-methods!)
% .plot        plotting flag:
%               true or false (boolean)
% 
% specific options for fft-based solvers: Subfields of structure OPTIONS:
% (required for fft-based solvers only)
% .tol         solver relative residual:
%               usually about 1e-10
% .maxit       solver maximum iteration number:
%               usually number of measurements or less
% .cond        solver regularization parameter:
%               usually 1e6, best performance in most cases according to experience
% .verbose     solver verbosity:
%               0, 1 or 2
% .kalstr      solver filter strength:
%               zero or some percent of variance
%               (approximate solution of kriging equations for no-zero value)
% .flag_Strang solver method flag:
%               true or false (boolean)
%               (true may perform better for highly correlated fields)
%
% copyright 2007 by Wolfgang Nowak.
% Affiliation:
%   Institute of Hydraulic Engineering (IWS)
%   Chair for Hydraulics and Modelling of Hydrosystems (LH2)
%   University of Stuttgart, Germany
% Email:
%   wolfgang.nowak@iws.uni-stuttgart.de
%
% Detailed information on the method used here:
% J. Fritz, W. Nowak and I. Neuweiler: FFT-based Algorithms for Kriging",
% submitted to Mathematical Geology (2007)
% 
% version 1.5,  01 august 2007 / WN

% general notation hints:
% s    unknowns
% y    observations
% b    base functions
% n    for number
% x1   spatial coordinate
% x2   spatial coordinate
% x3   spatial coordinate
% _len length
% d    grid spacing
% h    effective spacing
% _e   embedded
% _pts for points
% _    for vectors specifying values for each direction
%      vectors with same name but without _
%      are usually the product of the _ vector

%--------------------------------------------
% Checking input
%--------------------------------------------

% ensure consistency of grid types and solver methods...
if isequal(y.gridtype,'regular')   && isequal(options.solver,'fft');
  options.solver = 'fft-reg';
end
if isequal(y.gridtype,'irregular') && isequal(options.solver,'fft');
  options.solver = 'fft-irreg';
end
if isequal(y.gridtype,'irregular') && isequal(options.estvar,'speedy');
  warning('GENERAL_KRIGING:input','.options: ESTVAR option SPEEDY not allowed for Y.GRIDTYPE option IRREGULAR. Setting ESTVAR to NONE.')
  disp('Try again with ESTVAR options NONE, FULL or ONE-POINT.')
  options.estvar='none';
end
if isequal(options.estvar,'speedy') && isequal(options.superpos,'standard')
  warning('GENERAL_KRIGING:input','.options: ESTVAR options SPEEDY is useless without SUPERPOS option FFT. Setting ESTVAR to NONE.')
  disp('Try again with ESTVAR options NONE, FULL, HYBRID or ONE-POINT.')
  options.estvar='none';
end

% checking input parameter
s            = check_structure_s(s);

%--------------------------------------------
% Computing auxiliary quantities
%--------------------------------------------

% generate regular grid data and coordinates for unknowns
s            = ndgrid_setup(s);

% % include 1D case as simplified 2D case
% if numel(s.n_pts)==1
%   s.n_pts = [s.n_pts 1];
%   s.d_pts = [s.d_pts 1];
% end
% if isfield(y,'n_pts') && numel(y.n_pts)==1
%   y.n_pts = [y.n_pts 1];
% end

% generate regular embedded covariance function (required also for non-spectral methods!)
se           = find_embedding(s,options);
Qse1         = generate_covariance_embedded_first_row(s,se);
% compute fftn(Qse1) once for all...
if isequal(options.superpos,'fft')
  fftnQse1   = fftn(Qse1);
  fftnQse1   = real(fftnQse1);
else
  fftnQse1   = [];
end

% complete measurement information
y.npts       = numel(y.values);
if ~isfield(y,'indices') && isequal(y.gridtype,'regular');
  y.nd         = s.nd;
  y.d_pts      = s.d_pts.*y.d_ratio;
  y.d_tot      = y.d_pts.*y.n_pts;
  y.x_last     = y.x_first + y.d_ratio.*y.n_pts-1;
  y.x_mid      = y.x_first + y.d_ratio.*round(y.n_pts/2-1/2);
  % generate regular grid data and coordinates for measurements (required only for fft-solver on regular grid)
  aux          = reshape(1:s.npts,[s.n_pts 1]);          % regular grid of s
  switch y.nd
    case 1
      y.midindex = aux(y.x_mid(1));
      aux        = aux(y.x_first(1):y.d_ratio(1):y.x_last(1));
    case 2
      y.midindex = aux(y.x_mid(1),y.x_mid(2));
      aux        = aux(y.x_first(1):y.d_ratio(1):y.x_last(1),...
                       y.x_first(2):y.d_ratio(2):y.x_last(2));
    case 3
      y.midindex = aux(y.x_mid(1),y.x_mid(2),y.x_mid(3));
      aux        = aux(y.x_first(1):y.d_ratio(1):y.x_last(1),...
                       y.x_first(2):y.d_ratio(2):y.x_last(2),...
                       y.x_first(3):y.d_ratio(3):y.x_last(3));
  end
  y.indices    = aux(:);                             % measurement indices in field of unknowns
end

% start up embedding for fft-reg solver
if isequal(options.solver,'fft-reg')
  y.model    = s.model;
  y.variance = s.variance;
  y.lambda   = s.lambda;
  y.nugget   = s.nugget;
  y.micro    = s.micro;
  ye         = find_embedding(y,options);
  Qye1       = generate_covariance_embedded_first_row(y,ye);
end

% generate coordinates for measurements from indices of unknowns
y.x_pts     = cell(s.nd,1);
for i=1:s.nd
  y.x_pts{i} = s.x_pts{i}(y.indices);
end

if isequal(options.solver,'standard') || isequal(options.estvar,'standard')
  Qyy         = generate_covariance_full(s,y.x_pts);    % measurement error is added later
end

%--------------------------------------------
% Setting up blocks for Kriging system of equations
%--------------------------------------------

switch b.model
  case 'zero'
    b.HX     = [];                                      % no additional rows and columns
    b.invQbb = [];                                      % no additional rows and columns
    y.rhs    = y.values;                                % use plain measurement values
    b.rhs    = [];                                      % no additional rows and columns
    b.n      = 0;                                       % ensure that input is correct
  case 'known'
    b.HX     = [];                                      % no additional rows and columns
    b.invQbb = [];                                      % no additional rows and columns
    y.rhs    = y.values-b.beta_pri;                     % subtract known mean from measurement values
    b.rhs    = [];                                      % no additional rows and columns
    b.n      = 0;                                       % ensure that input is correct
  case 'uncertain'
    b.HX     = b.trends(y.indices);                     % base function values at measurement locations
    b.invQbb = inv(b.Qbb);                              % invQbb from Qbb
    y.rhs    = y.values;                                % use plain measurement values
    b.rhs    = -b.invQbb*b.beta_pri;                    % additional rows of prior mean coefficients
  case 'unknown'
    b.HX     = b.trends(y.indices);                     % base function values at measurement locations
    b.invQbb = zeros([b.n 1]);                          % invQbb is zero
    y.rhs    = y.values;                                % use plain measurement values
    b.rhs    = zeros([b.n 1]);                          % additional rows of zeros
  otherwise
    error('GENERAL_KRIGING.input.b: options for type of mean must be ZERO, KNOWN, UNKNOWN or UNCERTAIN')
end
R            = speye(y.npts)*y.error;                   % measurement error

%--------------------------------------------
% Solve Kriging system of equations
%--------------------------------------------

switch options.solver
  case 'standard'
    [aux1,aux2] = solve_kriging(Qyy ,R,b,y,s.n_pts,options);
  case 'fft-reg'
    [aux1,aux2] = solve_kriging(Qye1,R,b,y,s.n_pts,options);
  case 'fft-irreg'
    [aux1,aux2] = solve_kriging(Qse1,R,b,y,s.n_pts,options);
  otherwise
    error('SOLVE_KRIGING.input.options: solver option must be STANDARD or FFT-REG or FFT-IRREG')
end
switch b.model
  case {'uncertain' 'unknown'}
    Pbb       = - inv(b.HX'*aux2 + b.invQbb);
    beta      = - Pbb * (aux2'*y.rhs + b.invQbb*b.beta_pri);
    ksi       = aux1 - aux2*beta;
  otherwise
    Pbb       = [];
    beta      = [];
    ksi       = aux1;
end
    

%--------------------------------------------
% Evaluating estimator
%--------------------------------------------

% fluctuating part
estimate = superposition(Qse1,ksi,y.indices,s.n_pts,se.n_pts,s.nd,options,fftnQse1);

% base function contributions
switch b.model
  case 'zero'
    % add nothing
  case 'known'
    estimate = estimate + b.beta_pri;
  case 'uncertain'
    estimate = estimate + reshape(b.trends*beta,[s.n_pts 1]);
  case 'unknown'
    estimate = estimate + reshape(b.trends*beta,[s.n_pts 1]);
  otherwise
    error('EVALUATE_ESTIMATE.input.b: model must be ZERO or KNOWN or UNCERTAIN or UNKNOWN')
end

%--------------------------------------------
% Evaluating estimation variance
%--------------------------------------------

% initializing estimation variance to prior variance
est_var       = ones([s.n_pts 1])*s.variance;

% contributions of the mean value
for i=1:b.n
  Qsy_zi      = superposition(Qse1,aux2(:,i),y.indices,s.n_pts,se.n_pts,s.nd,options,fftnQse1);
  X_i         = reshape(b.trends(:,i),[s.n_pts 1]);
  for j=1:b.n
    if i~=j
      Qsy_zj  = superposition(Qse1,aux2(:,j),y.indices,s.n_pts,se.n_pts,s.nd,options,fftnQse1);
      X_i     = reshape(b.trends(:,j),[s.n_pts 1]);
    else
      Qsy_zj  = Qsy_zi;
      X_j     = X_i;
    end
    est_var     = est_var - (Qsy_zi-X_i).*Pbb(i,j).*(Qsy_zj-X_j);
  end
end

% dealing with measurement contributions
switch options.estvar
  case 'full'
    for i=1:y.npts
      y.rhs    = unit_vector(y.npts,i);
      switch options.solver
        case 'standard'
          aux1 = solve_kriging(Qyy ,R,b,y,s.n_pts,options,'aux1_only');
        case 'fft-reg'
          aux1 = solve_kriging(Qye1,R,b,y,s.n_pts,options,'aux1_only');
        case 'fft-irreg'
          aux1 = solve_kriging(Qse1,R,b,y,s.n_pts,options,'aux1_only');
      end
      Qsy_i    = shiftaround(Qse1,y.indices(i),s.n_pts,s.nd);
      S_i      = superposition(Qse1,aux1,y.indices,s.n_pts,se.n_pts,s.nd,options,fftnQse1);
      est_var  = est_var - S_i.*Qsy_i;
    end
  case 'one-point'
    Qse1_star       = Qse1.*Qse1;
    y.rhs           = ones([y.npts 1])*s.variance^2/(s.variance+y.error);
      switch options.solver
        case 'standard'
          Qyy_star  = Qyy.*Qyy;
          aux1      = solve_kriging(Qyy_star ,0,b,y,s.n_pts,options,'aux1_only');
        case 'fft-reg'
          Qye1_star = Qye1.*Qye1;
          aux1      = solve_kriging(Qye1_star,0,b,y,s.n_pts,options,'aux1_only');
        case 'fft-irreg'
          aux1      = solve_kriging(Qse1_star,0,b,y,s.n_pts,options,'aux1_only');
      end
    est_star        = superposition(Qse1_star,aux1,y.indices,s.n_pts,se.n_pts,s.nd,options,fftnQse1);
    est_var         = est_var - est_star;
  case 'speedy'
    y.mid_y  = find(y.indices == y.midindex);
    y.rhs    = unit_vector(y.npts,y.mid_y);
    switch options.solver
      case 'standard'
        aux1 = solve_kriging(Qyy ,R,b,y,s.n_pts,options,'aux1_only');
      case 'fft-reg'
        aux1 = solve_kriging(Qye1,R,b,y,s.n_pts,options,'aux1_only');
      case 'fft-irreg'
        aux1 = solve_kriging(Qse1,R,b,y,s.n_pts,options,'aux1_only');
    end
    [~,Se_i] = superposition(Qse1,aux1,y.indices,s.n_pts,se.n_pts,s.nd,options,fftnQse1);
    [~,Se1]  = shiftaround(Se_i,y.midindex,s.n_pts,s.nd,-1);
    Ze1      = Se1.*Qse1;
    clear X Se_i Se1 Qse1 aux1 aux2
    est_star = superposition(Ze1,ones([y.npts 1]),y.indices,s.n_pts,se.n_pts,s.nd,options,real(fftn(Ze1)));
    est_var  = est_var - est_star;
%   case 'hybrid'
%     warning('GENERAL_KRIGING:input','.options: ESTVAR options HYBRID is not implemented yet. Returning NaN.')
%     est_var  = zeros([s.n_pts 1])*NaN; % return all NaN to avoid misinterpretation
  case 'none'
    est_var  = zeros([s.n_pts 1])*NaN; % return all NaN to avoid misinterpretation
  otherwise
    error('GENERAL_KRIGING.input.options: estvar option must be FULL, ONE-POINT, SPEEY or HYBRID')
end

%-----------------------------------------
% plotting the results
%-----------------------------------------

if options.plot == true
  if ~(isempty(est_var)) && ~any(isnan(est_var(:)))
    subplot(2,1,1)
    plotter_nd(s.x_vec{1},s.x_vec{2},s.x_vec{3},estimate,'Kriging estimate',[1 1 1],s.n_pts,[],[],0,100,0);
    subplot(2,1,2)
    plotter_nd(s.x_vec{1},s.x_vec{2},s.x_vec{3},est_var,'Kriging variance',[1 1 1],s.n_pts,[],[],0,100,0,[1 0 0]);
  else
    plotter_nd(s.x_vec{1},s.x_vec{2},s.x_vec{3},estimate,'Kriging estimate',[1 1 1],s.n_pts,[],[],0,100,0);
  end
end






