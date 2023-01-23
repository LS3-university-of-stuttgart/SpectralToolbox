function [model,s,t]=randomize_model(model)
% RANDOMIZE_MODEL randomizes a geostatistical model MODEL
% version: 21 Nov 2009
% written by Wolfgang Nowak
% email: Wolfgang.Nowak@iws.uni-stuttgart.de
%
% RANDOMIZE_MODEL provides random structural parameters (normal or log-normal)
% with mean values provided in MODEL.THETA
% and covariances provided in MODEL.QTT
% and distribution type specified in MODEL.THETADIST
%
% required fields in input structure MODEL:
% .QTT       covariance matrix (symmetric, nonnegative)
% .N_THETA   number of uncertain parameters
% .THETA     mean values of parameters
% .THETADIST distribution type "pseudonormal" or "lognormal" or "none"
%
% output parameters:
% MODEL     randomized geostatistical model
% S         random Gaussian vector of parameter values
% .THETA    contains the tranfrormed (physical) values following .THETADIST
%           also returned in T

% draw random structural parameters for all uncertain ones
if any(model.Qtt~=0)
  switch model.thetadist
    case 'lognormal'
      % compute square root decomposition
      [V,D]       = eig(model.Qtt);
      sqrt_Qtt    = V*sqrt(D)*V';
      % correct between mean before and after exponentiation
      theta_mean  = model.theta - 0.5*diag(model.Qtt);
      % generate random values
      s           = sqrt_Qtt'*randn(model.n_theta,1);
      % take exponential to get desired distribution
      model.theta = exp(theta_mean + s);
    case 'pseudonormal'
      % in this mode, "s" AND "theta" are pseudo-normal
      % compute square root decomposition
      [V,D]       = eig(model.Qtt);
      sqrt_Qtt    = V*sqrt(D)*V';      
      % generate random values
      s           = sqrt_Qtt'*randn(model.n_theta,1);
      % produce correlated unit random numbers
      u           = normcdf_(s./sqrt(diag(model.Qtt)),0,1);
      % transform to bounded pseudo-normal variable
      CV               = sqrt(diag(model.Qtt))./model.theta;
      if CV > 2
        error('do not ask for CVs greater than 2 when expecting positive pseudo-normal numbers!')
      end
      betascale        = max(1+6*CV,2); % upper bound of beta distribution
      if model.n_theta > 0
        [p,ha,hb]      = betapdf_ms(0.2,1/betascale(1),1/betascale(1)*CV(1));
        model.theta(1) = model.theta(1)*betascale(1)*betainv(u(1),ha,hb);
      end
      if model.n_theta > 1
        [p,ha,hb]      = betapdf_ms(0.2,1/betascale(2),1/betascale(2)*CV(2));
        model.theta(2) = model.theta(2)*betascale(2)*betainv(u(2),ha,hb);
      end
      if model.n_theta > 2
        [p,ha,hb]      = betapdf_ms(0.2,1/betascale(3),1/betascale(3)*CV(3));
        model.theta(3) = model.theta(3)*betascale(3)*betainv(u(3),ha,hb);
      end
      if model.n_theta > 3
        [p,ha,hb]      = betapdf_ms(0.2,1/betascale(4),1/betascale(4)*CV(4));
        model.theta(4) = model.theta(4)*betascale(4)*betainv(u(4),ha,hb);
      end
    case 'none'
      % in this mode, the model remains deterministic
    otherwise
      error('randomize_model.m: model.thetadist: unknown type')
  end

  % store random theta values to corresponding fields of "model"
  if model.n_theta > 0
    model.variance = model.theta(1);
  end
  if model.n_theta > 1
    model.lambda   = model.theta(2:numel(model.lambda)+1)';
  end
  if model.n_theta > numel(model.lambda)+1
    model.kappa    = model.theta(numel(model.lambda)+2);
  end
  if model.n_theta > numel(model.lambda)+2
    error('randomize_model.m: model.n_theta too high')
  end

else
  % user did not ask for uncertainty in theta, so the model remains deterministic
  s = zeros(model.n_theta,1);
end

t = model.theta;

