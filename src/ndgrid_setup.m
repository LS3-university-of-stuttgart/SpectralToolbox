function s = ndgrid_setup(s)
% NDGRID_SETUP generates regular grid from basic grid information
% First node is at coordinate origin.
% Currently restricted to one to three dimensions in call to NDGRID_ND.
% 
% version 29 january 2008 / WN
%
% required input parameters:
% s.n_pts       vector of grid size per direction
% s.d_pts       vector of grid spacing per direction
%
% output parameters:
% s.x_pts       cell array of grid coordinates, each sized prod(nd)
% s.x_vec       cell array with grid vectors, each length n_pts
% s.nd          number of dimensions (scalar)
% s.d_tot       ND times one vector with total grid length

% checking required input parameters
if ~isfield(s,'n_pts')
  error('NDGRID_SETUP.input: input parameter S must contain a field S.N_PTS.')
end
if ~isfield(s,'d_pts')
  error('NDGRID_SETUP.input: input parameter S must contain a field S.D_PTS.')
end
if numel(s.n_pts) ~= numel(s.d_pts) || numel(s.n_pts) < 0 || numel(s.n_pts) > 3
  error('NDGRID_SETUP.input: size of s.n_pts and s.d_pts must both equal to the number of desired dimensions, i.e., 1,2,3.')
end

% accomplishing grid info fields
s.npts    = prod(s.n_pts);
s.nd      = numel(s.n_pts);
s.allpts  = reshape_nd(1:s.npts,s.n_pts);
s.d_tot   = s.d_pts .* s.n_pts; 

% defining domain coordinate vectors
s.x_vec   = cell(s.nd,1);
for i=1:s.nd
  if isfield(s,'n_add') % for the case of embedding
    s.x_vec{i}  = [(0:floor(s.n_pts(i)/2)), -ceil(s.n_pts(i)/2)+1:-1] *s.d_pts(i);
  else
    s.x_vec{i}  = (0:s.n_pts(i)-1)*s.d_pts(i);
  end
end
for i=s.nd+1:3
  s.x_vec{i}  = [];
end

% generating mesh of coordinates
s.x_pts   = ndgrid_nd(s.x_vec,s.nd);
