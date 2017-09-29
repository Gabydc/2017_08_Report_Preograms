
function [x, error, ii, flag] = pcg_m(A, x, b, M, max_it, tol)

%  -- Iterative template routine --
%     Univ. of Tennessee and Oak Ridge National Laboratory
%     October 1, 1993
%     Details of this algorithm are described in "Templates for the
%     Solution of Linear Systems: Building Blocks for Iterative
%     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
%     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
%     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
%
%  [x, error, iter, flag] = cg(A, x, b, M, max_it, tol)
%
% cg.m solves the symmetric positive definite linear system Ax=b 
% using the Conjugate Gradient method with preconditioning.
%
% input   A        REAL symmetric positive definite matrix
%         x        REAL initial guess vector
%         b        REAL right hand side vector
%         M        REAL preconditioner matrix
%         max_it   INTEGER maximum number of iterations
%         tol      REAL error tolerance
%
% output  x        REAL solution vector
%         error    REAL error norm
%         iter     INTEGER number of iterations performed
%         flag     INTEGER: 0 = solution found to tolerance
%                           1 = no convergence given max_it

  flag = 0;                                 % initialization
  ii = 0;
x0 = x;
  bnrm2 = norm( b );
  if  ( bnrm2 == 0.0 ), bnrm2 = 1.0; end
 
  r = b - A*x;
  error = norm( r ) / bnrm2;
  if ( error < tol ) return, end

  for ii = 1:max_it                       % begin iteration

     z  = M \ r;
     rho = (r'*z);

     if ( ii > 1 ),                       % direction vector
        beta = rho / rho_1;
        p = z + beta*p;
     else
        p = z;
     end

     q = A*p;
     alpha = rho / (p'*q );
     x = x + alpha * p;                    % update approximation vector

     r = r - alpha*q;                      % compute residual
     error = norm( r ) / bnrm2;            % check convergence
     if ( error <= tol ), break, end 

     rho_1 = rho;

     
     

        ee = norm(x - x0)/norm(x0);    
    trn = norm( b - A * x );
    error_r = norm(r)/norm(b);
     figure(120)
     color=[0.9 0.8 0.2];
     hline=plot(ii,log(norm(r)),'*','Color',color);
     hold on
       %  figure(120)
     color=[0.2 0.8 0.2];
     hline=plot(ii,log(trn),'s','Color',color);
     hold on
     color=[0.2 0.5 0.5];
     hline=plot(ii,log(error_r),'o','Color',color);
%      hold on
%      color=[0.8 0.3 0.5];
%      hline=plot(ii,log(),'d','Color',color);
%      hold on
%      
%     
     
     
  end
legend(['norm( r ) = ' num2str(norm(r))], ['norm( b - A * x ) =' ...
     num2str(trn)],['norm(r)/norm(b) = ' num2str(error_r)]);

  if ( error > tol ) flag = 1; end         % no convergence

% END cg.m