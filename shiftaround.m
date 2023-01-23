function [Qsy_i,Qsye_i] = shiftaround(Qse1,center,n_pts,nd,direction)
% SHIFTAROUND generates single lines of Qsy from first row of embedded Qss
% version 01 august 2007 / WN

if nargin < 5 || isempty(direction)
  direction = +1;
end

switch nd
  case 1
    [k1      ] =               center ;
    Qsye_i     = circshift(Qse1,direction*( k1       -1));
  case 2
    [k1,k2   ] = ind2sub(n_pts,center);
    Qsye_i     = circshift(Qse1,direction*([k1,k2   ]-1));
  case 3
    [k1,k2,k3] = ind2sub(n_pts,center);
    Qsye_i     = circshift(Qse1,direction*([k1,k2,k3]-1));
  otherwise
end
Qsy_i = extraction(Qsye_i,n_pts,nd);