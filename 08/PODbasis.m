%This program gives the matrix for the solution of the Laplace problem, the
%linear system is solved with the CGF conjugate gradient function

function [U,S]=PODbasis(x)
% Number of unknowns 
lx=size(x,1);
% Number of snapshots
ly=size(x,2);

% Normalize the snapshots
for i=1:ly
    x(:,i)=x(:,i)/norm(x(:,i)); 
end
%Compute mean of the snapshots
xm = mean(x,2);
% if ly == 2
% else
% for i=1:ly
%     x(:,i)=x(:,i)-xm;
% end
% end
D=(x'*x);
D = sparse(D);
[V,L] = eigs(D,ly);
S=zeros(lx,ly);
g=diag(L);
g=sqrt(abs(g));
S(1:ly,1:ly)=diag(g);
S=sparse(S);
%size(S)
V=sparse(V);
%size(V)
x=sparse(x);
%size(x)
U=x*V*S';
% size(U)

