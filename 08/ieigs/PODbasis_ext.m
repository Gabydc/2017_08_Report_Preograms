%This program gives the matrix for the solution of the Laplace problem, the
%linear system is solved with the CGF conjugate gradient function
clear all
%function [U,S]=PODbasis_ex(x,xn)
 x = [
     1    2   
     3    4   
     5    6   
     7    8   ];
 xn = [
     5      
     6       
     2       
     9  ];
 xn = xn/norm(xn);
 
 
% Number of unknowns 
lx=size(x,1);
% Number of snapshots
ly=size(x,2);

% % Normalize the snapshots
for i=1:ly
    x(:,i)=x(:,i)/norm(x(:,i));
    
end
x1 = [x xn];
%Compute mean of the snapshots
xm = mean(x,2);

xm1 = mean(x1,2)

if ly == 2
else
for i=1:ly
    x(:,i)=x(:,i)-xm
end
end




D=(1/ly)*x'*x;
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
U = x*V*S';
U = U(:,1:ly);
S = S(1:ly,:);
% size(U)
xax=1:ly;
figure
plot(log(g1),'ob')
%hold on
%plot(log(sqrt(abs(g))),'*r')
axis('tight')
%title(['Eigenvalues R=(1/m) X*X^T'],'FontSize',16);
ylabel('log(Value) ','FontSize',16)
            xlabel('Eigenvalue','FontSize',16)
            
            
% % Update the mean
 xmu = (1/(ly+1))*(ly*xm+xn);
% 
 gp = U'*(xn-xmu);
% 
 e = xn -(U*gp + xmu)
% 
 en = norm(e);
% 
if en == 0
    e = 0;
else
    e= e/en;
end

gamma = e'*(xn-xmu)

M0 = zeros(size(S))

 D1 = (ly/ly+1) * [S  gp*0] + (ly/(ly+1)^2) * [gp*gp'  gamma*gp]
 D2 = (ly/ly+1) * [gp'*0 0] + (ly/(ly+1)^2) * [gamma*gp' gamma^2]

D = [D1
    D2]
Uu = [U e];


        
            
            
            

