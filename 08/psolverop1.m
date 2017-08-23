

%p0 = x.pressure;

    switch lsolver
        case 1
            mrstModule add agmg
            solver = AGMGSolverAD('tolerance', tol);
        case 2
            solver = GMRES_ILUSolverAD('tolerance', tol);
        case 3
            solver = PCG_ICSolverAD_cn('tolerance', tol,'maxIterations', maxIterations,'x0',p0);
        case 4
            solver = DPCG_ICSolverAD_cn('tolerance', tol,'maxIterations', maxIterations, 'Z',Z,'cn',cn);
        otherwise
            solver = BackslashSolverAD();
    end