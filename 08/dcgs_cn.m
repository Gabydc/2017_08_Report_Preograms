
function [xf,flag,relres,ii,resvec,conda,conde] = dcgs_cn(A,b,Z,tol,maxit,M1,M2,x0,varargin)
%relres=1000000;
%DCGS   Deflated Conjugate Gradients Squared Method.
%   X = DCGS(A,B,Z) attempts to solve the system of linear equations A*X=B for
%   X, deflated. The N-by-N coefficient matrix A must be square and the right hand
%   side column vector B must have length N.
%
%   X = DCGS(AFUN,B,Z) accepts a function handle AFUN instead of the matrix A.
%   AFUN(X) accepts a vector input X and returns the matrix-vector product
%   A*X. In all of the following syntaxes, you can replace A by AFUN.
%
%   X = DCGS(A,B,Z,TOL) specifies the tolerance of the method. If TOL is []
%   then CGS uses the default, 1e-6.
%
%   X = CGS(A,B,Z,TOL,MAXIT) specifies the maximum number of iterations. If
%   MAXIT is [] then DCGS uses the default, min(N,20).
%
%   X = DCGS(A,B,Z,TOL,MAXIT,M) and X = DCGS(A,B,Z,TOL,MAXIT,M1,M2) use the
%   preconditioner M or M=M1*M2 and effectively solve the system
%   A*inv(M)*X = B for X. If M is [] then a preconditioner is not
%   applied.  M may be a function handle returning M\X.
%
%   X = DCGS(A,B,Z,TOL,MAXIT,M1,M2,X0) specifies the initial guess. If X0 is
%   [] then DCGS uses the default, an all zero vector.
%
%   [X,FLAG] = DCGS(A,B,...) also returns a convergence FLAG:
%    0 DCGS converged to the desired tolerance TOL within MAXIT iterations.
%    1 DCGS iterated MAXIT times but did not converge.
%    2 preconditioner M was ill-conditioned.
%    3 DCGS stagnated (two consecutive iterates were the same).
%    4 one of the scalar quantities calculated during DCGS became too
%      small or too large to continue computing.
%
%   [X,FLAG,RELRES] = DCGS(A,B,...) also returns the relative residual
%   NORM(B-A*X)/NORM(B). If FLAG is 0, then RELRES <= TOL.
%
%   [X,FLAG,RELRES,ITER] = DCGS(A,B,...) also returns the iteration number
%   at which X was computed: 0 <= ITER <= MAXIT.
%
%   [X,FLAG,RELRES,ITER,RESVEC] = DCGS(A,B,...) also returns a vector of the
%   residual norms at each iteration, including NORM(B-A*X0).
%
%   Example:
%      n = 21; A = gallery('wilk',n);  b = sum(A,2);
%      tol = 1e-12;  maxit = 15; M = diag([10:-1:1 1 1:10]);
%      x = cgs(A,b,tol,maxit,M);
%   Or, use this matrix-vector product function
%      %-----------------------------------------------------------------%
%      function y = afun(x,n)
%      y = [0; x(1:n-1)] + [((n-1)/2:-1:0)'; (1:(n-1)/2)'].*x+[x(2:n); 0];
%      %-----------------------------------------------------------------%
%   and this preconditioner backsolve function
%      %------------------------------------------%
%      function y = mfun(r,n)
%      y = r ./ [((n-1)/2:-1:1)'; 1; (1:(n-1)/2)'];
%      %------------------------------------------%
%   as inputs to CGS:
%      x1 = cgs(@(x)afun(x,n),b,tol,maxit,@(x)mfun(x,n));
%
%   Class support for inputs A,B,M1,M2,X0 and the output of AFUN:
%      float: double
%
%   See also BICG, BICGSTAB, BICGSTABL, GMRES, LSQR, MINRES, PCG, QMR,
%   SYMMLQ, TFQMR, ILU, FUNCTION_HANDLE.

%   Copyright 1984-2013 The MathWorks, Inc.

if (nargin < 3)
    error(message('MATLAB:dpcg:NotEnoughInputs'));
end

% Determine whether A is a matrix or a function.
[atype,afun,afcnstr] = iterchk(A);
if strcmp(atype,'matrix')
    % Check matrix and right hand side vector inputs have appropriate sizes
    [m,n] = size(A);
    xf = zeros(n,1);   
    if (m ~= n)
        error(message('MATLAB:dpcg:NonSquareMatrix'));
    end
    if ~isequal(size(b),[m,1])
        error(message('MATLAB:dpcg:RSHsizeMatchCoeffMatrix', m));
    end
else
    m = size(b,1);
    n = m;
    if ~iscolumn(b)
        error(message('MATLAB:dpcg:RSHnotColumn'));
    end
end

% Check matrix Z has appropriate size
[mz,nz] = size(Z);
if ~isequal(size(b),[mz,1])
    error(message('MATLAB:dpcg:WrongDeflationMAtrixSize', mz));
end

% Assign default values to unspecified parameters
if (nargin < 4) || isempty(tol)
    tol = 1e-6;
end
warned = 0;
if tol < eps
    warning(message('MATLAB:dpcg:tooSmallTolerance'));
    warned = 1;
    tol = eps;
elseif tol >= 1
    warning(message('MATLAB:dpcg:tooBigTolerance'));
    warned = 1;
    tol = 1-eps;
end
if (nargin < 5) || isempty(maxit)
    maxit = min(n,20);
end

% Check for all zero right hand side vector => all zero solution
n2b = norm(b);                     % Norm of rhs vector, b
if (n2b == 0)                      % if    rhs vector is all zeros
    xf = zeros(n,1);                % then  solution is all zeros
    flag = 0;                      % a valid solution has been obtained
    relres = 0;                    % the relative residual is actually 0/0
    ii = 0;                      % no iterations need be performed
    ee = 0;
    resvec = 0;                    % resvec(1) = norm(b-A*x) = norm(0)
    warning(message('MATLAB:dpcg: Solution is all zeros (rhs is zeros)'));
    if (nargout < 2)
        itermsg('dpcg',tol,maxit,0,flag,ii,NaN);
    end
    return
end


if ((nargin >= 6) && ~isempty(M1))
    existM1 = 1;
    [m1type,m1fun,m1fcnstr] = iterchk(M1);
    if strcmp(m1type,'matrix')
        if ~isequal(size(M1),[m,m])
            error(message('MATLAB:dpcg:WrongPrecondSize', m));
        end
    end
else
    existM1 = 0;
    m1type = 'matrix';
end
if ((nargin >= 7) && ~isempty(M2))
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

if ((nargin >= 8) && ~isempty(x0))
    if ~isequal(size(x0),[n,1])
        error(message('MATLAB:pcg:WrongInitGuessSize', n));
    else
        x = x0;
    end
else
    x = zeros(n,1);
end

if ((nargin > 8) && strcmp(atype,'matrix') && ...
        strcmp(m1type,'matrix') && strcmp(m2type,'matrix'))
    error(message('MATLAB:pcg:TooManyInputs'));
end

E = Z'*A*Z;
EI = sparse(inv(E));
if n<5000
 Q = Z*EI*Z';
Pd= eye(n) - A*Q;
  Pda=Pd*A;   
[Va,Da] = eigs(inv(M1)*Pda*inv(M2),n);
 Da=diag(Da);
 Da=abs(Da(nz+1:n));
 digits(4)
   conda=max(Da)/min(Da);
  [Ve,De] = eigs(E);
    De=diag(De);
  De=abs(De);
  digits(4)
  conde=max(De)/min(De);
 %conda=cond(inv(M1)*Pda*inv(M2),2);
%conde=condest(E);
 end
% Set up for the method
[pb]=defvec(Z,EI,A,b);
lb=M1\b;
plb=M1\pb;
nor=abs(lb'*lb);
r=b-iterapp('mtimes',afun,atype,afcnstr,x,varargin{:});
[r0]=defvec(Z,EI,A,r);
flag = 1;
xmin = x;                          % Iterate which has minimal residual so far
imin = 0;                          % Iteration at which xmin was computed
tolb = tol * n2b;                  % Relative tolerance
r = b - iterapp('mtimes',afun,atype,afcnstr,x,varargin{:});
normr = norm(r);                   % Norm of residual
normr_act = normr;

if (normr <= tolb)                 % Initial guess is a good enough solution
    flag = 0;
    relres = normr / n2b;
    ii = 0;
    resvec = normr;
    ee=0;
    xf=0;
     warning(message('MATLAB:dpcg: Initial guess is a good enough solution'));
    if (nargout < 2)
        itermsg('pcg',tol,maxit,0,flag,ii,relres);
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
for ii=1:maxit
    ee = 0;
    if ((rho == 0) || isinf(rho))
        flag = 4;
        warning(message('MATLAB:dpcg: residual is too large'));
        break
    end
    % Compute alpha=(r0,r0)/(PAp0,p0);
    q=iterapp('mtimes',afun,atype,afcnstr,p0,varargin{:});
    [ap]=defvec(Z,EI,A,q);
    pq = p0' * ap;
    
    
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
        y = iterapp('mldivide',m1fun,m1type,m1fcnstr,ap,varargin{:});
        if ~all(isfinite(r0))
            flag = 2;
            return
        end
    else % no preconditioner
        y = ap;
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
        warning(message('MATLAB:dpcg: residual is too large or too small'));
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
    pr=tdefvec(Z,EI,A,r);
    ee=abs(pr'*pr)/nor;
    %      figure(152)
    %      hl1=semilogy(iter,ee,'p','Color',color);
    %      hold on
    flag=0;
     x=xf;
 
    relres(ii)=ee;
    if (ee>=tol)
        flag=1;
        
    end
    if flag==0
      %  fprintf('DICCG Converged in %d iterations\n',ii);
        
        break
    end
   
end


[xf]=tdefvec(Z,EI,A,xf);
[qb]=qvec(Z,EI,b);
xf=qb+xf;
end
function[qx]=qvec(z,ei,x)
qx=z'*x;
qx=ei*qx;
qx=z*qx;
end
function[px]=defvec(z,ei,a,x)
[qx]=qvec(z,ei,x);
px=x-a*qx;
end
function[px]=tdefvec(z,ei,a,x)
ax=a'*x;
[qax]=qvec(z,ei,ax);
px=x-qax;
end



