function ue = padding(u,n_pts,n_pts_e,nd)
% PADDING does an embedding in larger field of zeros
% Currently restircted to 1, 2 or 3 dimensions.
% version 29 january 2008 / WN
%
% input parameters:
% u           field to be embedded with NDIM(UE) = number of dimensions
% n_pts_e     size of field,
%             must have same number of dimensions as U,
% n_pts_e     size of larger field,
%             must have same number of dimensions as UE,
%             must be larger or equal to SIZE(UE)
% nd          number of dimensions
%
% output parameters:
% u           extracted field

if nd ~= numel(n_pts_e) || nd < 0 || nd > 3
  error('PADDING.input: SIZE(U) and NDIM(N_PTS_E) must be equal and 1, 2, or 3.')
end
if any(n_pts_e < size(u))
  error('PADDING.input: N_PTS_E must be greater or equal to SIZE(U) in each element.')
end

ue    = zeros([n_pts_e,1]);
switch nd
  case 1
    ue(1:n_pts(1))                       = u;
  case 2
    ue(1:n_pts(1),1:n_pts(2))            = u;
  case 3
    ue(1:n_pts(1),1:n_pts(2),1:n_pts(3)) = u;
end