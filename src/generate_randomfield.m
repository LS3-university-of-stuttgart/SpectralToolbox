function [Y,beta,FFTQe] = generate_randomfield(s,b,options,white_noise,FFTQe)


% GENERATE_RANDOMFIELD generates realizations of random space variables.
% version: 2008-03-04
% written by Wolfgang Nowak
% email: Wolfgang.Nowak@iws.uni-stuttgart.de
%
% GENERATE_FIELD uses the spectral method (Dietrich & Newsam / Kitanidis)
% to generate multi-Gaussian correlated random fields.
% if FFTQe is recycled from previous calls, it speeds up
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
% b.model          : uncertain, known, zero or unknown
% b.beta_pri       : uncertain mean value (scalar). default value is zero.
% b.Qbb            : variance of uncertain mean value (scalar, set this to zero for known mean)
%                    default value is zero.
% s.micro          : microscale smoothing parameter relative to lambda
%                    (scalar, after Kitanidis book). default value is
%                    zero
% s.nugget         : additional nugget effect (salar, relative to model.variance)
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

% Reference for Dietrich & Newsam method:
% @article{Dietrich_Newsam_93,
%    author    = {C. R. Dietrich and G. N. Newsam},
%    title     = {A Fast and Exact Method For Multidimensional {Gaussian} Stochastic Simulations},
%    journal   = {Water Resour. Res.},
%    pages     = {2861-2869},
%    year      = {1993},
%    volume    = {29},
%    number    = {8}
% }
% @article{Dietrich_Newsam_97,
%    author    = {C. R. Dietrich and G. N. Newsam},
%    title     = {Fast and Exact Simulation of Stationary {Gaussian} Processes
%                 Through Circulant Embedding of the Covariance Matrix},
%    journal   = {SIAM J. ScI. Comput.},
%    pages     = {1088-1107},
%    year      = {1997},
%    volume    = {18},
%    number    = {4}
% }
%
% Reference for nice explanation of embedding methods:
% @article{Nowak_al_2003,
%    author    = {W. Nowak and S. Tenkleve and O. A. Cirpka},
%    title     = {Efficient computation of linearized cross-covariance and
%                 auto-covariance matrices of interdependent quantities},
%    journal   = {Math. Geol.},
%    pages     = {53-66},
%    year      = {2003},
%    volume    = {35},
%    number    = {1}
% }
%
% Reference for micro_D
% @book{Kitanidis_book,
%    author    = {P. K. Kitanidis},
%    title     = {Introduction to Geostatistics},
%    year      = {1997},
%    publisher = {Cambridge University Press},
%    address   = {Cambridge}
% }
% Reference for Kitanidis method: personal correspondence (sorry...)
%
% Reference for Zinn & Harvey-type fields
% @article{Zinn_Harvey_2003,
%    author    = {B. Zinn and C. F. Harvey},
%    title     = {When good statistical models of aquifer heterogeneity go bad:
% 	A comparison of flow, dispersion, and mass transfer in connected and multivariate
% 	{Gaussian} hydraulic conductivity fields},
%    journal   = {Water Resour. Res.},
%    year      = {2003},
%    volume    = {39},
%    number    = {3},
%    pages     = {1051, doi:10.1029/2001WR001146}
% }
% checking input parameter
if nargin == 0
  error('GENERATE_RANDOMFIELD:missing_input','input parameters (model and grid structures) are missing.')
elseif nargin > 5
  error('GENERATE_RANDOMFIELD:too_much_input','only two input parameters (model and grid structures) required.')
end

% checking required input parameters
s = check_structure_s(s);
s = ndgrid_setup(s);
b = check_structure_b(b);

% recycling FFTQe, if provided

if nargin ~=5 || isempty(FFTQe)
  [FFTQe,se] = initialize_FFT_cov(s,options);
else
  se         = find_embedding(s,options);
end


% computing eigenvalues of decompositions (due to sqrt)
s.sqrtQlambda  = sqrt(FFTQe/se.npts);

% loop to ensure non-NaN and non-Inf fields
while 1==1

  % computing embedded random field
  if s.flag_kit == false    % dietrich and newsam method
    epsilon        = complex(white_noise.a,white_noise.b);
    % epsilon(1)     = 0; % use only for zero mean of periodic fields, do not use for finite domains!
    Ye             = real(ifftn(epsilon.*s.sqrtQlambda))*se.npts;
  else                    % kitanidis method
    epsilon        = exp(1i*angle(fftn(white_noise.a)));
    epsilon(1)     = 0;   % this MUST be done because in the Kitanidis-type generation, abs(epsilon(1)) would otherwise always be unity, adding a bias to the mean!
    Ye             = real(ifftn(epsilon.*s.sqrtQlambda))*se.npts;
  end

  if s.flag_zh > 0    % producing connected fields (zinn and harvey)
    Ye             = Ye/sqrt(s.variance);  % ensuring unit variance
    % Do not correct this field to enforce zero-mean:
    % Ye must be zero mean in the ensemble, but enforcing zero mean
    % in single non-ergodic realizations messes up the mean after transformation!
    if s.flag_zh == 1 % low -K inclusions
      Ye           = -erfinv(2*(1-s.zh_smoother)*erf(abs(Ye)*sqrt(0.5))-(1-s.zh_smoother))*sqrt(2);
    elseif s.flag_zh == 2 % high-K inclusions
      Ye           = +erfinv(2*(1-s.zh_smoother)*erf(abs(Ye)*sqrt(0.5))-(1-s.zh_smoother))*sqrt(2);
    end
    Ye             = Ye*sqrt(s.variance);  % re-installing desired variance
  end

  % extracting random field from embedded field
  Y              = extraction(Ye,s.n_pts,s.nd);

  if s.flag_kit == 2
    % ensuring a correct normal distribution of mean values
    % !must be used for Kitanidis method and non-periodic fields!
    % otherwise, mean value would have a bimodal distribution (for some weird reasons which are too long to be explained here).
    % !! This correction is necessary, but it leads to spatially varying variance of Y in the ensemble mean !!
    % The correction is based on the relation (in matrix notation):
    % var(meanY) = E[u'*(Y-Yprime)*(Y-Yprime)'*u] = u'*CovYY*u;
    u              = ones(s.n_pts)/s.npts;
    Qsu            = convolution_FFT(s.FFTQse,u);
    uQsu           = u'*Qsu;
    Y              = Y - mean(Y(:)) + white_noise.c*sqrt(uQsu);
  end

  % adding (uncertain) mean value
  % If the field is too small to be ergodic, there is an additional uncertainty in the mean value
  % !! This effect is correct and desired !!
  switch b.model
    case 'zero'
      % nothing is necessary!
    case 'known'
      beta  = b.beta_pri;
    case 'uncertain'
      beta  = b.beta_pri + chol(b.Qbb)'*white_noise.d;
    case 'unknown'
      error('GENERATE_RANDOMFIELD.input.b: options UNKNOWN inly possible for kriging, but not for field generation')
    otherwise
      error('GENERATE_RANDOMFIELD.input.b: options for type of mean must be ZERO, KNOWN, UNKNOWN or UNCERTAIN')
  end
  Y       = Y + reshape(b.trends*beta,s.n_pts);

  % exit only when all non-NaN and all non-Inf
  if any(isnan(Y(:)))==false && any(isinf(Y(:)))==false, break, end
end