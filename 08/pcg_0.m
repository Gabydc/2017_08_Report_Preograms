function [x,flag,relres,ii,resvec,conda] = pcg_0(A,b,tol,maxit,M1,M2,x0,varargin)
r0 = b - A*x0;
rn0 = norm(r0);
r0 = M1\r0;
p0 = M2\r0;
figure(120)
clf
if norm(r0) ~= 0
    for ii = 1: maxit
        w = A * p0;
        alpha = (r0' * r0)/(p0' * w);
        x = x0 + alpha * p0;
        r = r0 - alpha * (M1 \ w);
    
        beta = (r' * r) / (r0' * r0);
        p = M2 \ r + beta * p0;
        flag = 0;
        rn = norm( r );
        resvec(ii) = rn;

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
        if (ee >= tol)
            flag = 1;
        end
        if flag == 0
            break
        end
        x0 = x;
        r0 = r;
        p0 = p;
    end
     legend(['norm( r ) = ' num2str(norm(r))], ['norm( b - A * x ) =' ...
     num2str(trn)],['norm(r)/norm(b) = ' num2str(error_r)]);
    relres = ee;
    conda = 0;
else
    x = x0;

    

 end