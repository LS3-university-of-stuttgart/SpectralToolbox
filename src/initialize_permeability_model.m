function [n_white_noise,sample_size,parameters] = initialize_permeability_model()

% defining domain discretization
s.n_pts        = [50 50];
s.d_pts        = [1   1  ];

% defining required parameters for geostatistical model
%s.model       = 'matern';               % model type
s.model       = 'gaussian';               % model type
s.variance    = 1;                      % variance for zero separation distance
s.lambda      = [20 20];                % correlation length

% defining optional parameters for geostatistical model
s.micro       = 0;                      % microscale smoothness parameter
s.nugget      = 0;                      % additional nugget effect (relative to variance)
s.kappa       = 0.5;                    % shape parameter for Matern covariance function
s.flag_kit    = 0;                      % true:use kitanidis method
s.flag_zh     = 1;                      % true: used connected zinn&harvey fields
s.zh_smoother = 0.01;                   % regularization parameter for zh fields
s.periodicity = [0 0];                  % conductivity fields not periodic in directions with 0
s.maxprime    = 7;                      % performance parameter for embedding size

% defining parameters for geostatistical trend model
b.model       = 'uncertain';            % uncertain, known, zero or unknown
b.n           = 1;                      % number of trend functions
b.beta_pri    = [log(1e-9)]';           % expected valuees of trend coefficients
b.Qbb         = 1;                      % (co-)variance of mean/trend coefficients

% this piece of code is dependent on model.nbeta!
% defining mean and trend functions
b.trend{1}    = ones(s.n_pts);          % constant mean function
b.trends      = b.trend{1}(:);

% specific options for fft-based routines
options.maxprime = 7;                   % embedding optimization parameter (2,3,5,7,...)

%% Randomize

%Computing in advance
se = find_embedding_and_check_inputs(s,options,b);
n_white_noise = compute_n_white_noise(se,s,b);

% Generate one Random field and Setup
disp([datestr(now) ' generating random field, including setup...'])
%Produce white noise
white_noise.all = randn(n_white_noise,1);
%Backtransformation of white noise
white_noise = reshape_white_noise(white_noise,s,b,se);
[Y,beta,FFTQe]    = generate_randomfield(s,b,options,white_noise);

sample_size = size(Y);

% For MCMC (not in use here)
parameters.s = s;
parameters.se = se;
parameters.b = b;
parameters.FFTQe = FFTQe;
parameters.options = options;

end
