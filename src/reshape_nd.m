function u = reshape_nd(u,n_pts)
% RESHAPE_ND evaluates matlab-function RESHAPE for all dimensionalities
% including the special case of 1D, based on a single size vector
% version 29 january 2008 / WN
% 
% input parameters:
% u          vector or array to be reshaped
% n_pts      size vector
%
% output parameters:
% u          reshaped result

if prod(n_pts) ~= numel(u)
  error('RESHAPE_ND.input: number of elements in U and size specified in N_PTS must be compatible.')
end

if numel(n_pts) == 1
  u = reshape(u,[n_pts 1]);
else
  u = reshape(u,n_pts);
end

