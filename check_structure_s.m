function s = check_structure_s(s)
% CHECK_STRUCTURE_S checks and initializes the variable S 
% version 29 january 2008 / WN
%
% Required geostatistical input parameters:
% s.model          : parametric geostatistical model name (string code)
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

% checking required input parameters
if ~isfield(s ,'d_pts')       || isempty(s.d_pts)
  error('GENERATE_Y:missing_input','s.d_pts not set.')
end
if ~isfield(s ,'n_pts')       || isempty(s.n_pts)
  error('GENERATE_Y:missing_input','s.n_pts not set.')
end
if ~isfield(s,'model')        || isempty(s.model)
  error('GENERATE_Y:missing_input','s.model not set.')
end
if ~isfield(s,'variance')     || isempty(s.variance)
  error('GENERATE_Y:missing_input','s.variance not set.')
end
if ~isfield(s,'lambda')       || isempty(s.lambda)
  error('GENERATE_Y:missing_input','s.lambda not set.')
end

% checking optional input parameters and setting defaults
if ~isfield(s ,'nd')          || isempty(s.nd)
  s.nd       = numel(s.n_pts);
end
if ~isfield(s ,'npts')        || isempty(s.npts)
  s.npts     = prod(s.n_pts);
end
if ~isfield(s,'micro')        || isempty(s.micro)
  s.micro    = 0;
end
if ~isfield(s,'nugget')       || isempty(s.nugget)
  s.nugget   = 0;
end
if ~isfield(s,'kappa')        || isempty(s.kappa)
  s.kappa    = 1.5;
end
if ~isfield(s,'periodicity')  || isempty(s.periodicity)
  s.periodicity = zeros(s.nd,1);
end
if ~isfield(s,'flag_kit')     || isempty(s.flag_kit)
  s.flag_kit = 0;
end
if ~isfield(s,'flag_zh')      || isempty(s.flag_zh)
  s.flag_zh  = 0;
end
if ~isfield(s,'zh_smoother')  || isempty(s.zh_smoother)
  s.zh_smoother = 0;
end

% checking correctness of some basic input properties
if ~isscalar(s.variance)      || s.variance < 0
  error('GENERATE_Y:variance:incorrect_input','variance must be a nonnegative scalar.')
end
if numel(s.lambda) ~= s.nd    || any(s.lambda) < 0
  error('GENERATE_Y:lambda:incorrect_input','lambda must be a nonnegative vector of same size as domain_len.')
end
if ~isscalar(s.micro)         || s.micro < 0
  error('GENERATE_Y:micro:incorrect_input','micro must be a nonnegative scalar.')
end
if numel(s.d_pts) ~= s.nd     || any(s.d_pts) <= 0
  error('GENERATE_Y:d_pts:incorrect_input','d_pts must be a nonnegative vector of same size as domain_len.')
end
if numel(s.n_pts) ~= s.nd     || any(s.n_pts) <= 0
  error('GENERATE_Y:n_pts:incorrect_input','n_pts must be a nonnegative vector of same size as domain_len.')
end
if numel(s.periodicity) ~= s.nd
  error('GENERATE_Y:periodicity:incorrect_input','periodicity must be a flag vector of same size as domain_len.')
end
if ~isscalar(s.flag_kit) || s.flag_kit ~= 0 && s.flag_kit ~= 1 && s.flag_kit ~= 2
  error('GENERATE_Y:flag_kit:incorrect_input','value of flag_kit must be 0, 1 or 2.')
end
if ~isscalar(s.flag_zh) || s.flag_zh ~= 0 && s.flag_zh ~= 1 && s.flag_zh ~= 2
  error('GENERATE_Y:flag_zh:incorrect_input','value of flag_zh must be 0, 1 or 2.')
end
if ~isscalar(s.zh_smoother) || s.zh_smoother < 0 || s.zh_smoother > 0.25
  error('GENERATE_Y:zh_smoother:incorrect_input','value of zh_smoother must be between 0 and 0.25.')
end
