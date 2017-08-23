%This program gives the matrix for the solution of the Laplace problem, the
%linear system is solved with the CGF conjugate gradient function

function [U,S]=PODbasis_ex(x)

% Normalize the snapshots
for i=1:size(x,2)
    x(:,i)=x(:,i)/norm(x(:,i));
end

%Compute mean of the snapshots
xm = mean(x,2);

% Number of unknowns 
lx=size(x,1);
% Number of snapshots
ly=size(x,2);


D=x'*x;
[V,L] = eigs(D,ly);
%size(D)
%size(V)
%size(L)
S=zeros(lx,ly);
g=diag(L);
g1=sqrt(abs(g));
%g1=1./g;
S(1:ly,1:ly)=diag(g1);
mg=max(g1);
S=sparse(S);
%size(S)
V=sparse(V);
%size(V)
x=sparse(x);
%size(x)
U=x*V*S';
% size(U)
xax=1:ly;
figure
plot(log(sqrt(abs(g1))),'ob')
hold on
plot(log(sqrt(abs(g))),'*r')
axis('tight')
title(['Eigenvalues R=X*X^T'],'FontSize',16);
ylabel('log(Value) ','FontSize',16)
            xlabel('Eigenvalue','FontSize',16)
            
            
% Update the mean
xmu = (1/(ly+1))*(m*xm+xn);

g = U'*(xn-xmu);

e = xn -(U*g + xmu);

en = norm(e);

if en == 0
    e = 0;
else
    e= e/en;
end


    
        
            
            
            

