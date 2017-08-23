function [x,flag,relres,iter,resvec,conda] = pcg_cn(A,b,tol,maxit,M1,M2,x0,varargin)
x=[];
flag=[];
relres=[];
iter=[];
resvec=[];
conda=[];
%PCG   Preconditioned Conjugate Gradients Method.
%   X = PCG(A,B) attempts to solve the system of linear equations A*X=B for
%   X. The N-by-N coefficient matrix A must be symmetric and positive
%   definite and the right hand side column vector B must have length N.
%
%   X = PCG(AFUN,B) accepts a function handle AFUN instead of the matrix A.
%   AFUN(X) accepts a vector input X and returns the matrix-vector product
%   A*X. In all of the following syntaxes, you can replace A by AFUN.
%
%   X = PCG(A,B,TOL) specifies the tolerance of the method. If TOL is []
%   then PCG uses the default, 1e-6.
%
%   X = PCG(A,B,TOL,MAXIT) specifies the maximum number of iterations. If
%   MAXIT is [] then PCG uses the default, min(N,20).
%
%   X = PCG(A,B,TOL,MAXIT,M) and X = PCG(A,B,TOL,MAXIT,M1,M2) use symmetric
%   positive definite preconditioner M or M=M1*M2 and effectively solve the
%   system inv(M)*A*X = inv(M)*B for X. If M is [] then a preconditioner
%   is not applied. M may be a function handle MFUN returning M\X.
%
%   X = PCG(A,B,TOL,MAXIT,M1,M2,X0) specifies the initial guess. If X0 is
%   [] then PCG uses the default, an all zero vector.
%
%   [X,FLAG] = PCG(A,B,...) also returns a convergence FLAG:
%    0 PCG converged to the desired tolerance TOL within MAXIT iterations
%    1 PCG iterated MAXIT times but did not converge.
%    2 preconditioner M was ill-conditioned.
%    3 PCG stagnated (two consecutive iterates were the same).
%    4 one of the scalar quantities calculated during PCG became too
%      small or too large to continue computing.
%
%   [X,FLAG,RELRES] = PCG(A,B,...) also returns the relative residual
%   NORM(B-A*X)/NORM(B). If FLAG is 0, then RELRES <= TOL.
%
%   [X,FLAG,RELRES,ITER] = PCG(A,B,...) also returns the iteration number
%   at which X was computed: 0 <= ITER <= MAXIT.
%
%   [X,FLAG,RELRES,ITER,RESVEC] = PCG(A,B,...) also returns a vector of the
%   estimated residual norms at each iteration including NORM(B-A*X0).
%
%   Example:
%      n1 = 21; A = gallery('moler',n1);  b1 = A*ones(n1,1);
%      tol = 1e-6;  maxit = 15;  M = diag([10:-1:1 1 1:10]);
%      [x1,flag1,rr1,iter1,rv1] = pcg(A,b1,tol,maxit,M);
%   Or use this parameterized matrix-vector product function:
%      afun = @(x,n)gallery('moler',n)*x;
%      n2 = 21; b2 = afun(ones(n2,1),n2);
%      [x2,flag2,rr2,iter2,rv2] = pcg(@(x)afun(x,n2),b2,tol,maxit,M);
%
%   Class support for inputs A,B,M1,M2,X0 and the output of AFUN:
%      float: double
%
%   See also BICG, BICGSTAB, BICGSTABL, CGS, GMRES, LSQR, MINRES, QMR,
%   SYMMLQ, TFQMR, ICHOL, FUNCTION_HANDLE.

%   Copyright 1984-2013 The MathWorks, Inc.

if (nargin < 2)
    error(message('MATLAB:pcg:NotEnoughInputs'));
end

% Determine whether A is a matrix or a function.
[atype,afun,afcnstr] = iterchk(A);
if strcmp(atype,'matrix')
    % Check matrix and right hand side vector inputs have appropriate sizes
    [m,n] = size(A);
    if (m ~= n)
        error(message('MATLAB:pcg:NonSquareMatrix'));
    end
    if ~isequal(size(b),[m,1])
        error(message('MATLAB:pcg:RSHsizeMatchCoeffMatrix', m));
    end
else
    m = size(b,1);
    n = m;
    if ~iscolumn(b)
        error(message('MATLAB:pcg:RSHnotColumn'));
    end
end

% Assign default values to unspecified parameters
if (nargin < 3) || isempty(tol)
    tol = 1e-6;
end
warned = 0;
if tol <= eps
    warning(message('MATLAB:pcg:tooSmallTolerance'));
    warned = 1;
    tol = eps;
elseif tol >= 1
    warning(message('MATLAB:pcg:tooBigTolerance'));
    warned = 1;
    tol = 1-eps;
end
if (nargin < 4) || isempty(maxit)
    maxit = min(n,20);
end

% Check for all zero right hand side vector => all zero solution
n2b = norm(b);                    % Norm of rhs vector, b
if (n2b == 0)                      % if    rhs vector is all zeros
    x = zeros(n,1);                % then  solution is all zeros
    flag = 0;                      % a valid solution has been obtained
    relres = 0;                    % the relative residual is actually 0/0
    iter = 0;                      % no iterations need be performed
    resvec = 0;                    % resvec(1) = norm(b-A*x) = norm(0)
    if (nargout < 2)
        itermsg('pcg',tol,maxit,0,flag,iter,NaN);
    end
    return
end

if ((nargin >= 5) && ~isempty(M1))
    existM1 = 1;
    [m1type,m1fun,m1fcnstr] = iterchk(M1);
    if strcmp(m1type,'matrix')
        if ~isequal(size(M1),[m,m])
            error(message('MATLAB:pcg:WrongPrecondSize', m));
        end
    end
else
    existM1 = 0;
    m1type = 'matrix';
end

if ((nargin >= 6) && ~isempty(M2))
    existM2 = 1;
    [m2type,m2fun,m2fcnstr] = iterchk(M2);
    if strcmp(m2type,'matrix')
        if ~isequal(size(M2),[m,m])
            error(message('MATLAB:pcg:WrongPrecondSize', m));
        end
    end
else
    existM2 = 0;
    m2type = 'matrix';
end

if ((nargin >= 7) && ~isempty(x0))
    if ~isequal(size(x0),[n,1])
        error(message('MATLAB:pcg:WrongInitGuessSize', n));
    else
        x = x0;
    end
else
    x = zeros(n,1);
end

if ((nargin > 7) && strcmp(atype,'matrix') && ...
        strcmp(m1type,'matrix') && strcmp(m2type,'matrix'))
    error(message('MATLAB:pcg:TooManyInputs'));                             
end

% Set up for the method
flag = 1;
xmin = x;                          % Iterate which has minimal residual so far
imin = 0;                          % Iteration at which xmin was computed
tolb = tol * n2b;                  % Relative tolerance
r0 = b - iterapp('mtimes',afun,atype,afcnstr,x,varargin{:});
normr = norm(r0);                   % Norm of residual
normr_act = normr;
lb=M1\b;
nor=abs(lb'*lb);

if n<5000

% [Va,Da] = eigs(inv(M1)*A*inv(M2),n);
%  Da=diag(Da);
%  Da=abs(Da);
%   conda=max(Da)/min(Da);
%  [Ve,De] = eigs(E,n);
%  De=diag(De);
%  De=abs(De);
%  conde=max(De)/min(De);
conda=condest(inv(M1)*A);
end

if (normr <= tolb)                 % Initial guess is a good enough solution
    flag = 0;
    relres = normr / n2b;
    iter = 0;
    resvec = normr;
    if (nargout < 2)
        itermsg('pcg',tol,maxit,0,flag,iter,relres);
    end
    return
end

resvec = zeros(maxit+1,1);         % Preallocate vector for norm of residuals
resvec(1,:) = normr;               % resvec(1) = norm(b-A*x0)
normrmin = normr;                  % Norm of minimum residual
rho = 1;
stag = 0;                          % stagnation of the method
moresteps = 0;
maxmsteps = min([floor(n/50),5,n-maxit]);
maxstagsteps = 3;


% Preconditioning the residual r0=M1\r0;
if existM1
    r0 = iterapp('mldivide',m1fun,m1type,m1fcnstr,r0,varargin{:});
    if ~all(isfinite(r0))
        flag = 2;
        return
    end
else % no preconditioner
    r0 = r0;
end
% Preconditioning the search direction p0=M2\r0;
if existM2
        p0 = iterapp('mldivide',m2fun,m2type,m2fcnstr,r0,varargin{:});
        if ~all(isfinite(p0))
            flag = 2;
            return
        end
    else % no preconditioner
        p0 = r0;
    end

%nor=1;
    rho1 = rho;
    rho = r0' * r0;
    
% loop over maxit iterations (unless convergence or failure)
for ii = 1 : maxit
     
    if ((rho == 0) || isinf(rho))
        flag = 4;
        break
    end
    % Compute alpha=(r0,r0)/(PAp0,p0);
    q=iterapp('mtimes',afun,atype,afcnstr,p0,varargin{:});
   
    pq = p0' * q;
  
    
    if ((pq <= 0) || isinf(pq))
        flag = 4;
        break
    else
        alpha = rho / pq;
    end
    if isinf(alpha)
        flag = 4;
        break
    end
     xf=x+alpha*p0;
     if existM1
    y = iterapp('mldivide',m1fun,m1type,m1fcnstr,q,varargin{:});
    if ~all(isfinite(r0))
        flag = 2;
        return
    end
else % no preconditioner
    y = q;
end
     r=r0-alpha*y;
     rho1 = rho;
    rho = r' * r;
    
%          if (ii == 1)
%         beta=(r'*r)/(r0'*r0);
%         p = p0;
%     else
        beta = rho / rho1;
        if ((beta == 0) || isinf(beta))
            flag = 4;
            break
        end
        p = M2\r + beta * p0;
%     end
    

% % Preconditioning the search direction p0=M2\r0;
% if existM2
%         p0 = iterapp('mldivide',m2fun,m2type,m2fcnstr,r,varargin{:});
%         if ~all(isfinite(p0))
%             flag = 2;
%             return
%         end
%     else % no preconditioner
%         p0 = r0;
%     end 

     

    % p=M2\r+beta*p0;  
     p0=p;
     r0=r;
     color=[0.1 0.5 0.5];
    normr = norm(r);
    normr_act = normr;
    resvec(ii+1,1) = normr;
      
      ee=abs(r'*r)/nor;
%      figure(152)
%      hl1=semilogy(iter,ee,'p','Color',color);
%      hold on
       flag=0;

     if (ee>=tol)
         flag=1;
     end
     if flag==0
         break
     end     
     x=xf;
    resvec(ii)=abs(r'*r); 
    relres(ii)=ee;
end

 iter=ii;
% figure(200)
% plot(x)
% hold on



% for ii = 1 : maxit
% 
% % returned solution is first with minimal residual
% if (flag == 0)
%     relres = normr_act / n2b;
% else
%     r_comp = b - iterapp('mtimes',afun,atype,afcnstr,xmin,varargin{:});
%     if norm(r_comp) <= normr_act
%         x = xmin;
%         iter = imin;
%         relres = norm(r_comp) / n2b;
%     else
%         iter = ii;
%         relres = normr_act / n2b;
%     end
% end
% 
% % truncate the zeros from resvec
% if ((flag <= 1) || (flag == 3))
%     resvec = resvec(1:ii+1,:);
% else
%     resvec = resvec(1:ii,:);
% end
% 
% % only display a message if the output flag is not used
% if (nargout < 2)
%     itermsg('pcg',tol,maxit,ii,flag,iter,relres);
% end
