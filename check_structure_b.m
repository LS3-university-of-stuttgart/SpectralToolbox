function b = check_structure_b(b)
% CHECK_STRUCTURE_B checks and initializes the variable B
% version 29 january 2008 / WN
%
% Required geostatistical input parameters:
% b.model          : uncertain, known, zero or unknown
% b.n              : number of basis functions
% b.beta_pri       : prior mean coefficients for trend functions
% b.Qbb            : uncertainty of prior mean: covariance matrix for trend coefficients
% b.trends         : cell array of trend functions (vector of ones for constant mean)
% b.trend          : array of trend functions


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
% s.flag_kit       : using kitanidis method to enforce spectrum.
%                    0 - dietrich and newsam method (default)
%                    1 - kitanidis method
%                    2 - kitanidis method with enforced gaussian distribution of mean value
% s.flag_zh:       : using zinn&harvey method to generate connectivity,
%                    0 - standard fields (default)
%                    1 - low- K inclusions
%                    2 - high-K inclusions
% s.zh_smoother:   : regularization parameter for zh fields.
%                    admissable range is between 0 (no effect) and 0.25


if ~isfield(b,'model')        || isempty(b.model)
  b.beta      = 'known';
end
if ~isfield(b,'n')            || isempty(b.n)
  b.n         = 1;
end
if ~isfield(b,'beta_pri')     || isempty(b.beta_pri)
  b.beta_pri  = 0;
end
if ~isfield(b,'Qbb')          || isempty(b.Qbb)
  b.Qbb       = 0;
end
if ~isfield(b,'trends')       || isempty(b.trends)
  b.trends{1} = 1;
end
if ~isfield(b,'trend')        || isempty(b.trend)
  b.trend     = 1;
end


% checking correctness of some basic input properties
if det(b.Qbb) < 0
  error('GENERATE_Y:Qbb:incorrect_input','Qbb must be nonnegative.')
end
if size(b.Qbb,1) ~= size(b.Qbb,2) || size(b.Qbb,1) ~= numel(b.beta_pri) || b.n ~= numel(b.beta_pri)
  error('GENERATE_Y:nbeta:incorrect_input','sizes and dimensions of Qbb, beta and nbeta must be compatible.')
end
if numel(b.trend) ~= b.n || b.n ~= round(b.n) || b.n < 0
  error('GENERATE_Y:n:incorrect_input','n.b bust be positive integer, and as large as trends cell array.')
end
