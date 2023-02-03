function [e] = find_embedding(s,options)
% FIND_EMBEDDING chooses minimum embedding size given grid and geostatistical model
% version 01 august 2007 / WN
%
% required input parameters:
% model_Y      : parametric geostatistical model
% variance_Y   : variance of model
% lambda_Y     : vector of correlation length in y,x,z-direction
% micro_Y      : microscale smoothing parameter relative to lambda_Y
%                (scalar, after Kitanidis book)
% d_tot        : length of domain in y,x,z-direction
% n_el         : number of elements in y,x,z-direction

% finding appropriate size for the embedding (approximately
% ensuring positive definiteness of embedded covariance matrix)
switch s.model
  case 'gaussian'
    e.minsize  = 3 + s.micro;
  case 'exponential'
    e.minsize  = 5 + s.micro;
  case 'spherical'
    e.minsize  = 1 + s.micro;
  case 'matern'
    e.minsize  = 5 + s.micro;
  otherwise
    error('GENERAL_KRIGING.find_embedding.input: geostatistical model must be GAUSSIAN or EXPONENTIAL or SPHERICAL')
end
e.maxsize = s.d_tot./s.lambda;

% defining embedded domain size
e.d_pts        = s.d_pts;
e.n_pts        = ceil(max(2*s.d_tot./s.d_pts , s.d_tot./s.d_pts + e.minsize.*s.lambda./s.d_pts));
e.n_pts        = nicer_primes(e.n_pts,options.maxprime);
e.npts         = prod(e.n_pts);
e.d_tot        = e.n_pts.*s.d_pts;
e.n_add        = e.n_pts - s.n_pts;
e.d_tot        = e.d_tot - s.d_tot;
e.x_pts        = cell(s.nd,1);

% generate full coordinate grid
e              = ndgrid_setup(e);

% ...periodicity of distances
% for i=1:s.nd
%   e.x_pts{i}   = min(e.x_pts{i},e.d_tot(i)-e.x_pts{i});
% end
% e.x_pts{1}((e.n_pts(1)/2+1):end,:) = -e.x_pts{1}((e.n_pts(1)/2+1):end,:);
% e.x_pts{2}(:,(e.n_pts(2)/2+1):end) = -e.x_pts{2}(:,(e.n_pts(2)/2+1):end);

% % smoothing the transition at the point of reflection using a microscale smoothing approach
% for i=1:s.nd
%   e.x_pts{i}  = sqrt((e.d_tot(i)/2).^2+0.1*s.lambda(i).^2)  -  sqrt((e.d_tot(i)/2-e.x_pts{i}).^2+0.1*s.lambda(i).^2);
%   e.x_pts{i}  = sqrt((e.d_tot(i)/2).^2+s.lambda(i).^2)  -  sqrt((e.d_tot(i)/2-e.x_pts{i}).^2+s.lambda(i).^2);
% end