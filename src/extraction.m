function u = extraction(ue,n_pts,nd)
% EXTRACTION reads smaller result from embedded field
% version 01 august 2007 / WN

switch nd
  case 1
    u = ue(1:n_pts(1));
  case 2
    u = ue(1:n_pts(1),1:n_pts(2));
  case 3
    u = ue(1:n_pts(1),1:n_pts(2),1:n_pts(3));
end