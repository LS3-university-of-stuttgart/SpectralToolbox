function ui = injection(u,H,n_pts)
% INJECTION writes values to specified locations into a field of zeros
% Like the POKE command in old BASIC
% version 29 january 2008 / WN
%
% input parameters:
% u          vector of values to be written into UI
% H          vector of indices (single-index, not subscripts)
%            to place values in U into field UI
% n_pts      size vector of desired field UI
%
% output parameters:
% ui         field of zeros, sized N_PTS, with values in U
%            written to locations H

if numel(u) ~= numel(H)
  error('INJECTION.input: NUMEL(U) and NUMEL(H) must be equal.')
end
if prod(n_pts) < max(H(:))
  error('INJECTION.input: indices in H must not exceed PROD(N_PTS).')
end

ui    = zeros(n_pts);
ui(H) = u;