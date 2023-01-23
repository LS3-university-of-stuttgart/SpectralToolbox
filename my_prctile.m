function y = my_prctile(x,p)
% Y = PRCTILE(X,P) returns the P-percentile of the values in the column vector X
% without drawing on statistics toolbox (may not be available on all systems)
% version 16 august 2006 / WN

x      = x(:);
n      = size(x,1);
if n==1
  y    = x;
else
  x      = sort(x,1);
  q      = [0 100*(0.5:(n-0.5))./n 100]';
  xx     = x([1,1:end,end]);
  y      = interp1q(q,xx,p);
end