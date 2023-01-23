function nicer = nicer_primes(number,pmax)
% NICER_PRIMES  looks for a larger number with prime factors
% smaller than pmax. This is useful for optimizing the embedding
% size in spectral methods to make sure the FFT is fast.
% The number is not necessarily the smallest one,
% but tends to be quite close...
%
% NICER = NICER_PRIMES(NUMBER) returns a number NICER >= NUMBER 
% which has nice prime factors. NUMBER must be a scalar, vector
% or matrix with positive integers. 
%
% NICER = NICER_PRIMES(NUMBER,PMAX) specifies the larges allowed
% prime factor in NICER. PMAX must be a single prime number.
% If no pmax is specified, NICER_PRIMES uses the default, PMAX = 7.
%
% Jochen Fritz, 2006.
% Email: jochen.fritz@iws.uni-stuttgart.de
% Version 1.2 by Wolfgang Nowak, 29 June 2007
% Email: wolfgang.nowak@iws.uni-stuttgart.de
%
% New in version 1.1:
% control of input parameters
%
% New in version 1.2:
% restructured code for better human readibility

if nargin < 1
  error('Not enough input arguments!')
end
if nargin > 2
  error('Too many input arguments!')
end
if nargin == 1 || isempty(pmax)
  pmax = 7;  
end
if any(number ~= round(number)) || any(number) <= 0
  error('number must be positive integers')
end
if numel(pmax) ~=1 || pmax < 2 || numel(factor(pmax)) ~= 1
  error('pmax must be a single prime number')
end

nicer           = ones(size(number));
for i = 1:numel(number),
  factors       = factor(number(i));
  while 1==1,
    nicefactors   = factors(factors <= pmax);
    bad_factors   = factors(factors >  pmax);
    nicer(i)    = prod(nicefactors) * nicer(i);
    rest        = prod(bad_factors);
    if rest == 1, break, end
    rest        = rest + 1;
    factors     = factor(rest);
  end
end

% search for possible smaller number based on smaller primes...
if pmax > 2
  pmax2  = max(primes(pmax-1)); % choose next smaller prime below pmax
  nicer2 = nicer_primes(number,pmax2); % check out value
  nicer  = min(nicer,nicer2);
end