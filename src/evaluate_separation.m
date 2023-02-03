function h_eff = evaluate_separation(x_pts,lambda,rotation,micro,nd)
% EVALUATE_SEPARATION evaluates effective separation distance given coordinates differences and scales
% version 01 august 2007 / WN

% evaluating effective separation distance
switch nd
  case 1
    h_eff      = x_pts{1}/lambda(1);
  case 2
    x_rotation =  x_pts{1}*cosd(rotation(1)) + x_pts{2}*sind(rotation(1));
    y_rotation =  -x_pts{1}*sind(rotation(1)) + x_pts{2}*cosd(rotation(1));
    h_eff      = sqrt((x_rotation/lambda(1)).^2 + (y_rotation/lambda(2)).^2);
  case 3
    x_rotation = x_pts{1}*(cosd(rotation(1))*cosd(rotation(3)) - cosd(rotation(2))*sind(rotation(1))*sind(rotation(3))) + ...
                 x_pts{2}*(sind(rotation(1))*cosd(rotation(3)) + cosd(rotation(2))*cosd(rotation(1))*sind(rotation(3))) + ...
                 x_pts{3}*sind(rotation(2))*sind(rotation(3));
    y_rotation = x_pts{1}*(-cosd(rotation(1))*sind(rotation(3)) - cosd(rotation(2))*sind(rotation(1))*cosd(rotation(3))) + ...
                 x_pts{2}*(-sind(rotation(1))*sind(rotation(3)) + cosd(rotation(2))*cosd(rotation(1))*cosd(rotation(3))) + ...
                 x_pts{3}*sind(rotation(2))*cosd(rotation(3));
    z_rotation = x_pts{1}*sind(rotation(2))*sind(rotation(1)) - x_pts{2}*sind(rotation(2))*cosd(rotation(1)) + x_pts{3}*cosd(rotation(2));
    h_eff      = sqrt((x_rotation/lambda(1)).^2 + (y_rotation/lambda(2)).^2 + (z_rotation/lambda(3)).^2);
end

% applying microscale smoothing
h_eff        = sqrt(h_eff.^2 + micro.^2)-micro;
