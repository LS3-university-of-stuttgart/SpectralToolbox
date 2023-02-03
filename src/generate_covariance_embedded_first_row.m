function Qse1 = generate_covariance_embedded_first_row(s,e)
% GENERATE_COVARIACE_EMBEDDED_FIRST_ROW generates first row of embedded covariance matrix
% version 01 august 2007 / WN

% generating coordinate differences between all points and first point
dx_pts = cell(s.nd,1);
for i=1:s.nd
  dx_pts{i} = e.x_pts{i}-e.x_pts{i}(1);
end

h_eff   = evaluate_separation(dx_pts,s.lambda,s.rotation,s.micro,s.nd);
Qse1    = evaluate_covariance(s,h_eff);