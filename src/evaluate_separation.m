function [h_eff] = evaluate_separation(x_pts,lambda,micro,nd)
% EVALUATE_SEPARATION evaluates effective separation distance given coordinates differences and scales
% version 01 august 2007 / WN

% evaluating effective separation distance
switch nd
  case 1
    h_eff      = x_pts{1}/lambda(1);
  case 2
    h_eff      = sqrt((x_pts{1}/lambda(1)).^2 + (x_pts{2}/lambda(2)).^2);
  case 3
    h_eff      = sqrt((x_pts{1}/lambda(1)).^2 + (x_pts{2}/lambda(2)).^2 + (x_pts{3}/lambda(3)).^2);
end

% applying microscale smoothing
h_eff        = sqrt(h_eff.^2 + micro.^2)-micro;