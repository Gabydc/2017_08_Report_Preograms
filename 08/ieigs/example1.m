close all
clear all
x = 0 : 1/24 : 1;
x1 = x-0.5;
t = 0 : 1/49 : 2;
t1 = t-1;
for i =1 : 25
    for j = 1 : 50
        z(i,j) = exp(-x(i)*t(j))+sin(x(i)*t(j));
    end
end
mesh(z)
 X =[
     1    2
     3    4
     5    6
     7    8];
%  for i=1:size(X,2)
%      x(:,i)=x(:,i)/norm(x(:,i));
%  end

 D=x'*x;
[U,S]=PODbasis_ext(X);
break
[U1,S1]=PODbasissvd(X);