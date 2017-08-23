

p0 = rSol.pressure;
n=nx*ny*nz;
for i=1:numel(W)
p0(n+i) = 0;
end
    switch lsolver
        case 1
            mrstModule add agmg
            solver = AGMGSolverAD('tolerance', tol);
        case 2
            solver = GMRES_ILUSolverAD('tolerance', tol);
        case 3
            solver = PCG_ICSolverAD_cn('tolerance', tol,'maxIterations', maxIterations,'cn',cn,'x0',p0);
        case 4
            solver = DPCG_ICSolverAD_cn('tolerance', tol,'maxIterations', maxIterations, 'Z',Z,'cn',cn,'x0',p0);
        otherwise
            solver = BackslashSolverAD();
    end