function [Qyy] = generate_covariance_full(s,x_pts)
% GENERATE_COVARIANCE_FULL generates full covariance matrix (measurements)
% version 01 august 2007 / WN

% get problem size
npts  = numel(x_pts{1});

% generating coordinate differences between all points
dx_pts = cell(s.nd,1);
for i=1:s.nd
  dx_pts{i} = x_pts{i}*ones(1,npts)-ones(npts,1)*x_pts{i}';
end

% evaluate separation distances and covariance matrix
h_eff = evaluate_separation(dx_pts,s.lambda,s.rotation,s.micro,s.nd);
Qyy   = evaluate_covariance(s,h_eff);