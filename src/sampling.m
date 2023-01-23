function u = sampling(ui,H)
% SAMPLING reads values from specified locations within a vector
% version 29 january 2008 / WN
%
% input parameters:
% ui      vector (or array) of numers to be samlped from
% H       list of indices where to sample from UI
%
% output parameters:
% u       list of sampled values


if max(H(:)) > numel(ui)
  error('SAMPLING.input: indices in H must not exceed NUMEL(UI).')
end

u(H) = ui;