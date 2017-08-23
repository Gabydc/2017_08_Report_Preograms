function [solver] = psolop(lsolver,rSol, tol, maxIterations,cn,varargin)
   opt = struct( 'wells', [],'Z',[]);
 opt = merge_options(opt, varargin{:});
pinit = rSol.pressure;
N = length(pinit);
W = opt.wells;
for w = 1:numel(W)
pinit(N+w) = 0;
end
Z = opt.Z;

    switch lsolver
        case 1
            mrstModule add agmg
            solver = AGMGSolverAD('tolerance', tol);
        case 2
            solver = GMRES_ILUSolverAD('tolerance', tol);
        case 3
            solver = PCG_ICSolverAD_cn('tolerance', tol,'maxIterations', maxIterations,'cn',cn,'x0',pinit);
        case 4
            solver = DPCG_ICSolverAD_cn('tolerance', tol,'maxIterations', maxIterations, 'Z',Z,'cn',cn,'x0',pinit);
        otherwise
            solver = BackslashSolverAD();
    end
end