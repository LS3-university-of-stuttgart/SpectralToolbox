function [x, relres, iter, flag, resvec, condition] = gsts_pcg(Ae, b, n_, H, tol, maxit, cond1, verbose, kalstr, x, flag_Strang)
%GSTS_PCG  GENERALIZED SPECTRAL TOEPLITZ SOLVER USING PRECONDITIONED CONJUGATE GRADIENT METHOD
%   X = GSTS_PCG(AE,B,N_,H) solves the subset A(H,H)*X(H)=B(H) of a symmetric (block) Toeplitz
%   system A*X=B for X(H) without fully assembling A(H) or even A.
%   N_ is a vector that specifies the size of a regular grid (in all space directions) that
%   includes all points listed in H. H is an index list that enumerates node indices on the
%   regular grid where the valuex of B(H) and X(H) are to be injected. The nodes in H
%   do not have to form a regular subgrid. Ae is the first row/column of a circulant matrix obtained
%   from circulant embedding of A in a periodic domain that is larger than N_.
%
%   For 1-D problems AE is a vector, for 2-D problems it is a matrix,
%   for 3-D problems, it is a 3rd order tensor. B must have the same number of
%   dimensions as Ae.
%
%   GSTS_PCG(AE,B,N_,H,TOL) specifies the tolerance of the method.
%   If TOL is [] then GSTS_PCG uses the default, 1e-6.
%
%   GSTS_PCG(AE,B,N_,H,TOL,MAXIT) specifies the maximum number of iterations.  
%   If MAXIT is [] then GSTS_PCG uses the default, 200.
%
%   GSTS_PCG(AE,B,N_,H,TOL,MAXIT,COND1) specifies the minimum condition. If the condition 
%   of the circulant matrix generated from AE is larger than COND1, then the 
%   main diagonal in the circulant matrix is amplified to meet COND1. COND1 must
%   be a scalar value larger than 1.
%   If COND1 is [] then GSTS_PCG uses the default, 1000.
%
%   GSTS_PCG(AE,B,HN_,,TOL,MAXIT,COND1,VERBOSE) specifies runtime output. If VERBOSE is 1,
%   then only a summary is displayed. If VERBOSE is 2, every single iteration step 
%   generates output. If VERBOSE is 3, there is no output. 
%   If VERBOSE is [] then GSTS_PCG uses the default, 1.
%
%   GSTS_PCG(AE,B,N_,H,TOL,MAXIT,COND1,VERBOSE,KALSTR) specifies strength of the Kalman filter
%   applied to X in which KALSTR is the white noise level removed from X. KALSTR = 0 means 
%   no filter. Positive non-zero values of KALSTR result in a filtered approximation of X 
%   and leads to faster and more stable convergence. KALSTR must be a scalar value larger
%   than zero.
%   If KALSTR is [] then GSTS_PCG uses the default, 0.
%
%   GSTS_PCG(AE,B,N_,H,TOL,MAXIT,COND1,VERBOSE,KALSTR,X) specifies the initial guess. Initial 
%   guess X must be of same size as B. 
%   If X is [] then GSTS_PCG uses the default, an all zero vector/matrix.
%
%   GSTS_PCG(AE,B,N_,H,TOL,MAXIT,COND1,VERBOSE,KALSTR,X,FLAG_STRANG) specifies whether to use
%   an embedded preconditioner (0) or the original-size Strang preconditioner.
%   The default is 0 for te embedded preconditioner. The Strang
%   preconditioner performes better for large correlation length scales.
%
%   [X,RELRES] = GSTS_PCG(AE,B,N_,H,TOL,MAXIT,COND1,VERBOSE,KALSTR,X) also returns the
%   relative residual NORM(Be-Ae*Xe)/NORM(B). 
%
%   [X,RELRES,ITER] = GSTS_PCG(AE,B,N_,H,TOL,MAXIT,COND1,VERBOSE,KALSTR,X) also returns 
%   the iteration number at which X was computed: 0 <= ITER <= MAXIT.
%
%   [X,RELRES,ITER,FLAG] = GSTS_PCG(AE,B,N_,H,TOL,MAXIT,COND1,VERBOSE,KALSTR,X) also 
%   returns a convergence FLAG:
%     1 trial solution fulfills TOL. 
%     0 GSTS_PCG converged to the desired tolerance TOL within MAXIT iterations.
%    -1 GSTS_PCG iterated MAXIT times but did not converge.
%    -2 RELRES = NaN.
%
%   [X,FLAG,RELRES,ITER,RESVEC] = GSTS_PCG(AE,B,N_,H,TOL,MAXIT,COND1,VERBOSE,KALSTR,X) 
%   also returns the error vector of the final iterate.
%
%   [X,FLAG,RELRES,ITER,RESVEC,CONDITION] = GSTS_PCG(AE,B,N_,H,TOL,MAXIT,COND1,VERBOSE,KALSTR,X)
%   also returns the condition of AE.
%
%   Sascha Tenkleve, 2003.
%   Email: s.tenkleve@debitel.net
%   $Revision: 1.1 $ $Date: 2006/11/30 17:33:25 $ by Wolfgang Nowak
%
%   Version 1.0:
%   based on sts_pcg version 1.7, then added injection scheme
%   for irregular grids and renamed to gsts_pcg
%
%   New in version 1.11:
%   now can use the Strang preconditioner



% Initial checks
if nargin < 4,      error('Not enough input arguments!'); end
if nargin > 11,     error('Too many input arguments!');   end
if  isempty(Ae),    error('Ae is empty!');
 elseif isnan(Ae),  error('Ae is NaN!');
 elseif ischar(Ae), error('Ae is char!');   
end  
if  isempty(b),     error('B is empty!');
 elseif isnan(b),   error('B is NaN!');
 elseif ischar(b),  error('B is char!');   
end  

n_e                 = size(Ae);
if length(n_) ~= length(n_e), error('embedding dimension mismatch!'); end
if length(n_) >=1 & ((n_(1) == 1)  & (n_e(1) ~= 1)), error('embedding of singleton first dimension unpermissible.'); end
if length(n_) >=2 & ((n_(2) == 1)  & (n_e(2) ~= 1)), error('embedding of singleton second dimension unpermissible'); end
if length(n_) >=3 & ((n_(3) == 1)  & (n_e(3) ~= 1)), error('embedding of singleton third dimension unpermissible'); end
if any(n_>n_e), error('Some embedded sizes are less than the non-embedded sizes!'); end
if prod(n_)<max(H(:)), error('injection list H'); end

% Assign default values to unspecified parameters
if nargin < 5  | isempty(tol),     tol     = 1e-6;           end
if nargin < 6  | isempty(maxit),   maxit   = 200;            end
if nargin < 7  | isempty(cond1),   cond1   = 1000;           end
if nargin < 8  | isempty(verbose), verbose = 1;              end
if nargin < 9  | isempty(kalstr),  kalstr  = 0;              end
if nargin < 10 | isempty(x),       x       = zeros(size(b)); end
if nargin < 11 | isempty(flag_Strang), flag_Strang = 0;      end

% Check parameters
if ((nargin >= 5) & ~isempty(tol))
   if ( ischar(tol)     | isnan(tol)     | tol     <= 0  | size(tol)   ~=1 ) 
       error('TOL must be scalar value larger than zero!');    
   end
end
if ((nargin >= 6) & ~isempty(maxit))
   if ( ischar(maxit)   | isnan(maxit)   | maxit   <  1  | size(maxit) ~=1 ) 
       error('MAXIT must be scalar value larger or equal 1!'); 
   end
end
if ((nargin >= 7) & ~isempty(cond1))
   if ( ischar(cond1)   | isnan(cond1)   | cond1   <= 1  | size(cond1) ~=1 ) 
       error('COND1 must be scalar value larger than 1!');     
   end
end
if ((nargin >= 8) & ~isempty(verbose))
   if ( ischar(verbose) | isnan(verbose) | verbose <  0  | verbose > 2 | size(verbose) ~=1 ) 
       error('VERBOSE must be scalar of value  0,1 or 2!');        
   end
end
if ((nargin >= 9) & ~isempty(kalstr))
   if ( ischar(kalstr)  | isnan(kalstr)  | kalstr  <  0  | size(kalstr) ~=1 )                
       error('KALSTR must be scalar value larger or equal zero!'); 
   end
end
if ((nargin ==10) & ~isempty(x))
   if ( ~isequal(size(x), size(b))       | ischar(x)     | isnan(x)  )                       
       error('Initial guess XO must be a non-NaN vector of same size as B!');        
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
resvec                  = b-Iconvolution(fftae,x,H,n_,n_e);
relres                  = sqrt(sum(resvec(:).*resvec(:)))/bno;
dx                      = Iconvolution(invfftRStrang,resvec,H,n_,n_e);
res                     = sum(resvec(:).*dx(:));
if relres <= tol
  flag = 1;
end

% core of algorithm
while 1==1
  % checking break criteria
  if relres <= tol,                  break; end %converged
  if isnan(relres) == 1,  flag = -2; break; end %relres = NaN
  if iter   >= maxit,     flag = -1; break; end %max iter exceeded
  iter                = iter + 1; 

  % updating trial solution
  qx                  = Iconvolution(fftae,dx,H,n_,n_e);
  alpha               = res/sum(dx(:).*qx(:));
  x                   = x + alpha * dx;
  
  % computing error vector
  if mod(iter,20)==0;
    resvec            = b-Iconvolution(fftae,x,H,n_,n_e);
  else
    resvec            = resvec - alpha * qx;
  end
  relres              = sqrt(sum(resvec(:).*resvec(:)))/bno;
  
  % conjugate gradient search direction
  sx                  = Iconvolution(invfftRStrang,resvec,H,n_,n_e);
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

function Qu = Iconvolution(fftae,u,H,n_,n_e);
% Iconvolution computes Q*u with circulant matrix from first row
% using an injection technique for irregular grids
% version 15 march 2006 / WN

% injection
ui    = zeros([n_ ,1]);
ui(H) = u;

% embedding
ue    = zeros([n_e,1]);
if length(n_)==1
  ue(1:n_(1))                 = ui;
elseif length(n_)==2
  ue(1:n_(1),1:n_(2))         = ui;
elseif length(n_)==3
  ue(1:n_(1),1:n_(2),1:n_(3)) = ui;
end

% convolution: Que   = real(ifftn(fftn(ue).*fftae));
% taken to pieces to save memory...
Que   = fftn(ue);
Que   = Que.*fftae;
Que   = ifftn(Que);
Que   = real(Que);

% extraction
if length(n_)==1
  Qui = Que(1:n_(1));
elseif length(n_)==2
  Qui = Que(1:n_(1),1:n_(2));
elseif length(n_)==3
  Qui = Que(1:n_(1),1:n_(2),1:n_(3));
end

% sampling
Qu    = Qui(H);

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