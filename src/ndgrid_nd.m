function [grid_cell] = ndgrid_nd(vec_cell,nd)
% NDGRID_ND evaluates matlab-function NDGRID for all dimensionalities 1-3
% version 29 january 2008 / WN
% 
% input parameters:
% vec        cell array of coordinate vectors
% nd         number of dimensions
%
% output parameters:
% grid       cell array containing the grid coordinates
%            for each dimension

if iscell(vec_cell) == false
  error('NDGRID_ND.input.vec: VEC_CELL must be a cell array with as many cells as dimensions.')
end
if nd < 1 || nd > 3
  error('NDGRID_ND.input.nd: ND or number of dimensions must be 1, 2 or 3.')
end

grid_cell = cell(nd,1);
switch nd
  case 1
    grid_cell{1}  = vec_cell{1}';
    grid_cell{2}  = [];
    grid_cell{3}  = [];
  case 2
    [grid_cell{1},grid_cell{2}]  = ndgrid(vec_cell{1},vec_cell{2});
    grid_cell{3}  = [];
  case 3
    [grid_cell{1},grid_cell{2},grid_cell{3}]  = ndgrid(vec_cell{1},vec_cell{2},vec_cell{3});
end