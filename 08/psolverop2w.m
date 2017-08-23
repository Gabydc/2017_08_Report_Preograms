

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
fn = @(A, b) solver.solveLinearSystem(A, b);
psolve = @(state) incompTPFA_g_o(state, G, hT, fluid, 'wells', W, 'MatrixOutput',true,'LinSolve', fn);
            
    %% Define transport solver
switch tsolver
    case 1
        % Explicit tranport solver
        tsolve  = @(state, dT, fluid) explicitTransport(state, G, dT, rock, fluid,  'wells', W,  'verbose', false);
    case 2
        % Implicit transport solver: try with one time step
        tsolve  = @(state, dT, fluid) implicitTransport(state, G, dT, rock, fluid,  'wells', W, 'Verbose', false);
end