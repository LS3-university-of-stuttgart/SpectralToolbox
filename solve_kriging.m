function [aux1,aux2] = solve_kriging(Q,R,b,y,n_se,options,flag_aux1_only)
% SOLVE_KRIGING solves Kriging system of equations using user-specified solver
% version 01 august 2007 / WN

if nargin < 7 || isempty(flag_aux1_only)
  flag_aux1_only = 'both';
else
  flag_aux1_only = 'aux1_only';
end  

switch options.solver
  case 'standard'
    if ~isequal(flag_aux1_only,'aux1_only') && (isequal(b.model,'unknown') || isequal(b.model,'uncertain'))
      aux1      = (Q+R)\y.rhs;
      aux2      = (Q+R)\b.HX;
    elseif isequal(b.model,'known') || isequal(b.model,'zero') || isequal(flag_aux1_only,'aux1_only')
      aux1      = (Q+R)\y.rhs;
      aux2      = [];
    end
  case 'fft-reg'
    Q(1)     = Q(1) + R(1);
    if ~isequal(flag_aux1_only,'aux1_only') && (isequal(b.model,'unknown') || isequal(b.model,'uncertain'))
      aux1 = sts_pcg(Q, reshape(y.rhs,[y.n_pts,1]), options.tol, options.maxit, options.cond, options.verbose, options.kalstr, [], options.flag_Strang);
      aux2 = sts_pcg(Q, reshape(b.HX ,[y.n_pts,1]), options.tol, options.maxit, options.cond, options.verbose, options.kalstr, [], options.flag_Strang);
      aux2 = aux2(:);
      aux1 = aux1(:);
    elseif isequal(b.model,'known') || isequal(b.model,'zero') || isequal(flag_aux1_only,'aux1_only')
      aux1 = sts_pcg(Q, reshape(y.rhs,[y.n_pts,1]), options.tol, options.maxit, options.cond, options.verbose, options.kalstr, [], options.flag_Strang);
      aux2 = [];
      aux1 = aux1(:);
    end
  case 'fft-irreg'
    Q(1)     = Q(1) + R(1);
    if numel(n_se) == 1, n_se = [n_se 1]; end
    if ~isequal(flag_aux1_only,'aux1_only') && (isequal(b.model,'unknown') || isequal(b.model,'uncertain'))
      aux1 = gsts_pcg(Q, reshape(y.rhs,y.npts,1), n_se, y.indices, options.tol, options.maxit, options.cond, options.verbose, options.kalstr, []);
      aux2 = gsts_pcg(Q, reshape(b.HX ,y.npts,1), n_se, y.indices, options.tol, options.maxit, options.cond, options.verbose, options.kalstr, []);
      aux2 = aux2(:);
      aux1 = aux1(:);
    elseif isequal(b.model,'known') || isequal(b.model,'zero') || isequal(flag_aux1_only,'aux1_only')
      aux1 = gsts_pcg(Q, reshape(y.rhs,y.npts,1), n_se, y.indices, options.tol, options.maxit, options.cond, options.verbose, options.kalstr, []);
      aux2 = [];
      aux1 = aux1(:);
    end
  otherwise
    error('SOLVE_KRIGING.input.options: solver option must be STANDARD or FFT-REG or FFT-IRREG')
end