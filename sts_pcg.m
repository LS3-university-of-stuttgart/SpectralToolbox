function [x, relres, iter, flag, resvec, condition] = sts_pcg(Ae, b, tol, maxit, cond1, verbose, kalstr, x, flag_Strang)
%STS_PCG   SPECTRAL TOEPLITZ SOLVER USING PRECONDITIONED CONJUGATE GRADIENT METHOD
%   X = STS_PCG(AE,B) solves the Toeplitz system of linear equations A*X=B
%   for X without fully assembling A. AE is the generator for a circulant 
%   matrix that is obtained from periodic embedding of the generator for A. 
%   For 1-D problems AE is a vector, for 2-D problems it is a matrix,
%   for 3-D problems, it is a 3rd order tensor. B must have the same number of
%   dimensions as Ae.
%
%   STS_PCG(AE,B,TOL) specifies the tolerance of the method.
%   If TOL is [] then STS uses the default, 1e-6.
%
%   STS_PCG(AE,B,TOL,MAXIT) specifies the maximum number of iterations.  
%   If MAXIT is [] then STS uses the default, 200.
%
%   STS_PCG(AE,B,TOL,MAXIT,COND1) specifies the minimum condition. If the condition 
%   of the circulant matrix generated from AE is larger than COND1, then the 
%   main diagonal in the circulant matrix is amplified to meet COND1. COND1 must
%   be a scalar value larger than 1.
%   If COND1 is [] then STS uses the default, 1000.
%
%   STS_PCG(AE,B,TOL,MAXIT,COND1,VERBOSE) specifies runtime output. If VERBOSE is 1,
%   then only a summary is displayed. If VERBOSE is 2, every single iteration step 
%   generates output. If VERBOSE is 0, there is no output. 
%   If VERBOSE is [] then STS uses the default, 0.
%
%   STS_PCG(AE,B,TOL,MAXIT,COND1,VERBOSE,KALSTR) specifies strength of the Kalman filter
%   applied to X in which KALSTR is the white noise level removed from X. KALSTR = 0 means 
%   no filter. Positive non-zero values of KALSTR result in a filtered approximation of X 
%   and leads to faster and more stable convergence. KALSTR must be a scalar value larger
%   than zero.
%   If KALSTR is [] then STS uses the default, 0.
%
%   STS_PCG(AE,B,TOL,MAXIT,COND1,VERBOSE,KALSTR,X) specifies the initial guess. Initial 
%   guess X must be of same size as B. 
%   If X is [] then STS_PCG uses the default, an all zero vector/matrix.
%
%   STS_PCG(AE,B,TOL,MAXIT,COND1,VERBOSE,KALSTR,X,FLAG_STRANG) specifies whether to use
%   an embedded preconditioner (0) or the original-size Strang preconditioner.
%   The default is 0 for te embedded preconditioner. The Strang
%   preconditioner performes better for large correlation length scales.
%
%   [X,RELRES] = STS_PCG(AE,B,TOL,MAXIT,COND1,VERBOSE,KALSTR,X,FLAG_STRANG) also returns the
%   relative residual NORM(Be-Ae*Xe)/NORM(B). 
%
%   [X,RELRES,ITER] = STS_PCG(AE,B,TOL,MAXIT,COND1,VERBOSE,KALSTR,X,FLAG_STRANG) also returns 
%   the iteration number at which X was computed: 0 <= ITER <= MAXIT.
%
%   [X,RELRES,ITER,FLAG] = STS_PCG(AE,B,TOL,MAXIT,COND1,VERBOSE,KALSTR,X,FLAG_STRANG) also 
%   returns a convergence FLAG:
%     1 trial solution fulfills TOL. 
%     0 STS_PCG converged to the desired tolerance TOL within MAXIT iterations.
%    -1 STS_PCG iterated MAXIT times but did not converge.
%    -2 RELRES = NaN.
%
%   [X,FLAG,RELRES,ITER,RESVEC] = STS_PCG(AE,B,TOL,MAXIT,COND1,VERBOSE,KALSTR,X,FLAG_STRANG) 
%   also returns the error vector of the final iterate.
%
%   [X,FLAG,RELRES,ITER,RESVEC,CONDITION] = STS_PCG(AE,B,TOL,MAXIT,COND1,VERBOSE,KALSTR,X,FLAG_STRANG)
%   also returns the condition of AE.
%
%   Sascha Tenkleve, 2003, Wolfgang Nowak, 2007
%   Email: wolfgang.nowak@iws.uni-stuttgart.de
%   $Revision: 1.19 $ $Date: 2006/11/31 17:33:25 $ by Wolfgang Nowak
%
%   New in version 1.15:
%   stabilized preconditioned steepest descent replaced by
%   preconditioned conjugate gradient method
%
%   New in version 1.16:
%   introduced subfunctions for embedding, extraction,
%   convolution and deconvolution
%
%   New in version 1.17:
%   all PCG operations performed on non-embedded level
%
%   New in version 1.18:
%   now can use the Strang preconditioner
%
%   New in version 1.19:
%   added quicker function return for trivial cases

% Initial checks
if nargin < 2,      error('Not enough input arguments!'); end
if nargin > 9,      error('Too many input arguments!');   end
if  isempty(Ae),    error('Ae is empty!');
 elseif isnan(Ae),  error('Ae is NaN!');
 elseif ischar(Ae), error('Ae is char!');   
end  
if  isempty(b),     error('B is empty!');
 elseif isnan(b),   error('B is NaN!');
 elseif ischar(b),  error('B is char!');   
end  

if isscalar(b)
  x         = b/Ae(1);
  relres    = 0;
  iter      = 0;
  flag      = 0;
  resvec    = 0;
  condition = 1;
  return
end

n_e                 = size(Ae);
n_                  = size(b);
if length(n_) ~= length(n_e), error('embedding dimension mismatch!'); end
% if length(n_) >=1 & ((n_(1) == 1)  & (n_e(1) ~= 1)), warning('embedding of singleton first dimension may be erroneous'); end
% if length(n_) >=2 & ((n_(2) == 1)  & (n_e(2) ~= 1)), warning('embedding of singleton second dimension may be erroneous'); end
% if length(n_) >=3 & ((n_(3) == 1)  & (n_e(3) ~= 1)), warning('embedding of singleton third dimension may be erroneous'); end
if any(n_>n_e), error('Some embedded sizes are less than the non-embedded sizes!'); end

% Assign default values to unspecified parameters
if nargin < 3 | isempty(tol),         tol         = 1e-6;           end
if nargin < 4 | isempty(maxit),       maxit       = 200;            end
if nargin < 5 | isempty(cond1),       cond1       = 1000;           end
if nargin < 6 | isempty(verbose),     verbose     = 0;              end
if nargin < 7 | isempty(kalstr),      kalstr      = 0;              end
if nargin < 8 | isempty(x),           x           = zeros(size(b)); end
if nargin < 9 | isempty(flag_Strang), flag_Strang = 0;              end

% Check parameters
if ((nargin >= 3) & ~isempty(tol))
   if ( ischar(tol)     | isnan(tol)     | tol     <= 0  | size(tol)   ~=1 ) 
       error('TOL must be scalar value larger than zero!');    
   end
end
if ((nargin >= 4) & ~isempty(maxit))
   if ( ischar(maxit)   | isnan(maxit)   | maxit   <  1  | size(maxit) ~=1 ) 
       error('MAXIT must be scalar value larger or equal 1!'); 
   end
end
if ((nargin >= 5) & ~isempty(cond1))
   if ( ischar(cond1)   | isnan(cond1)   | cond1   <= 1  | size(cond1) ~=1 ) 
       error('COND1 must be scalar value larger than 1!');     
   end
end
if ((nargin >= 6) & ~isempty(verbose))
   if ( ischar(verbose) | isnan(verbose) | verbose <  0  | verbose > 2 | size(verbose) ~=1 ) 
       error('VERBOSE must be scalar of value  0,1 or 2!');        
   end
end
if ((nargin >= 7) & ~isempty(kalstr))
   if ( ischar(kalstr)  | isnan(kalstr)  | kalstr  <  0  | size(kalstr) ~=1 )                
       error('KALSTR must be scalar value larger or equal zero!'); 
   end
end
if ((nargin >= 8) & ~isempty(x))
   if ( ~isequal(size(x), [n_])          | ischar(x)     | isnan(x) )                       
       error('Initial guess X must be a non-NaN vector of same size as B!');        
   end
end
if ((nargin == 9) & ~isempty(x))
   if ( numel(flag_Strang)~=1)           | (flag_Strang ~= 0 & flag_Strang ~= 1)
       error('FLAG_STRANG must be boolean.');        
   end
end

% preparing Strang preconditioner
if flag_Strang == 1
  Strang                  = extraction(Ae,n_+1,n_e);
  for i=1:ndims(b);
    Strang                = Strang+flipdim(Strang,i);
  end
  Strang                  = extraction(Strang,n_,n_+1);
end

% preparing fftn(Ae) and fftn(Strang) for later recycling
fftw('planner','patient')
fftae                   = complex(real(fftn(Ae    ) + kalstr));
if flag_Strang == 1
  fftStrang             = complex(real(fftn(Strang) + kalstr));
else
  fftStrang             = fftae;
end
fftw('planner','estimate')

% regularization and inversion of Strang preconditioner
lmin                    = min(real(fftStrang(:)));
lmax                    = max(real(fftStrang(:)));
delta                   = max((lmax-cond1*lmin)/(cond1-1),eps);
invfftRStrang           = 1./(fftStrang + delta);
if lmin ~= 0
  condition = lmax/lmin;
else
  condition = inf;
end

% initialize algorithm
flag                    = 0;
iter                    = 0; 
bno                     = sqrt(sum(b(:).*b(:)));
if bno == 0,
  x(:)   = 0;
  relres = 0;
  flag   = 1; %RHS is all zero!
else
  resvec                = b-convolution(fftae,x);
  relres                = sqrt(sum(resvec(:).*resvec(:)))/bno;
  if relres <= tol
    flag = 1; %initial guess is good
  else
  dx                    = convolution(invfftRStrang,resvec);
  res                   = sum(resvec(:).*dx(:));
  end 
end

% core of algorithm
while 1==1
  % checking break criteria
  if flag   == +1,        flag = +1; break; end %no need to even start!
  if relres <= tol,       flag = +0; break; end %converged
  if isnan(relres) == 1,  flag = -2; break; end %relres = NaN
  if iter   >= maxit,     flag = -1; break; end %max iter exceeded
  iter                = iter + 1; 

  % updating trial solution
  qx                  = convolution(fftae,dx);
  alpha               = res/sum(dx(:).*qx(:));
  x                   = x + alpha * dx;
  
  % computing error vector
  if mod(iter,20)==0;
    resvec            = b-convolution(fftae,x);
  else
    resvec            = resvec - alpha * qx;
  end
  relres              = sqrt(sum(resvec(:).*resvec(:)))/bno;
  
  % conjugate gradient search direction
  sx                  = convolution(invfftRStrang,resvec);
  resold              = res;
  res                 = sum(sx(:).*resvec(:));
  beta                = res/resold;
  dx                  = sx + beta*dx;
  
  if (verbose == 2)
    fprintf('relres: %5.5e \t iter: %d \t alpha: %5.5e\n', relres, iter, alpha);
  end
end % end of core

if (verbose == 1)
    fprintf('relres: %5.5e \t iter: %d \t alpha: %5.5e\n', relres, iter, alpha);
end

% ----------------------------------------------------------
%   S U B F U N C T I O N   1
% ----------------------------------------------------------

function Qu = convolution(fftae,u);
% convolution computes Q*u with circulant matrix from first row
% version 15 march 2006 / WN

n_  = size(u);
n_e = size(fftae);

% embedding
ue  = zeros([n_e,1]);
if length(n_)==1
  ue(1:n_(1))                 = u;
elseif length(n_)==2
  ue(1:n_(1),1:n_(2))         = u;
elseif length(n_)==3
  ue(1:n_(1),1:n_(2),1:n_(3)) = u;
end

% convolution: Que   = real(ifftn(fftn(ue).*fftae));
% taken to pieces to save memory...
Que   = fftn(ue);
Que   = Que.*fftae;
Que   = ifftn(Que);
Que   = real(Que);

% extraction
if length(n_)==1
  Qu = Que(1:n_(1));
elseif length(n_)==2
  Qu = Que(1:n_(1),1:n_(2));
elseif length(n_)==3
  Qu = Que(1:n_(1),1:n_(2),1:n_(3));
end

% ----------------------------------------------------------
%   S U B F U N C T I O N   2
% ----------------------------------------------------------
function u = extraction(ue,n_,n_e);
% extracts finite domain from periodic domain
% version 15 march 2006 / WN

if length(n_)==1
  u = ue(1:n_(1));
elseif length(n_)==2
  u = ue(1:n_(1),1:n_(2));
elseif length(n_)==3
  u = ue(1:n_(1),1:n_(2),1:n_(3));
end