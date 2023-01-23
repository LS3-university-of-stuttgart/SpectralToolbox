function Q = evaluate_covariance(s,h_eff)
% EVALUATE_COVARIANCE evaluates covariance function for a list of separation distances
% version 01 august 2007 / WN

switch s.model
  case 'gaussian'
    Q        = s.variance * exp(-h_eff.^2);
  case 'exponential'
    Q        = s.variance * exp(-h_eff);
  case 'spherical'
    Q        = s.variance * (1 - 1.5*h_eff + 0.5*h_eff.^3);
    Q(h_eff>1) = 0;
  case 'matern'
    Q        = s.variance/(2^(s.kappa-1)*gamma(s.kappa))...
               *  (h_eff*sqrt(s.kappa)*2).^(s.kappa)...
               .*  besselk(s.kappa,h_eff*sqrt(s.kappa)*2);
    Q(h_eff==0) = s.variance;
  otherwise
    error('GENERAL_KRIGING.evaluate_covariance.input: geostatistical model must be GAUSSIAN or EXPONENTIAL or SPHERICAL')
end

% adding nugget effect for zero separation distance
if isfield(s,'nugget') && s.nugget ~= 0
  Q(h_eff==0) = Q(h_eff==0) + s.nugget;
end