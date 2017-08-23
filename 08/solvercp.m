close all
clear all
sz=5;
nx = sz; ny = sz; nz = 1;
nw=0;

G         = cartGrid([nx ny nz],[nx ny nz]*meter);
G         = computeGeometry(G);
rock.perm = repmat(100*milli*darcy, [G.cells.num, 1]);
rock.poro = repmat(0.3            , [G.cells.num, 1]);
per=0;
 lsize=round(sz*sz/sz); 
for i=1:2:sz
 rock.perm(1+lsize*(i-1):lsize*i)  = repmat(10^(-per), [lsize, 1])*100*milli*darcy();
end

for i=1:sz

    Z((1:sz)+sz*(i-1),i)=1;
end


x = linspace(0, 1, 11) .';
y = linspace(1, 0, 11) .';
figure(100)
plotCellData(G, rock.perm);
hold on
%We define the relative permeability and the capillary pressure in form of tables, and let the relative permeability curves be quadratic and the capillary function linear. The strength of the capillary pressure is decided by cap_scale. The capillary pressure is defined in the non-wetting phase, p_c = p_nw - p_w.
hT = simpleComputeTrans(G, rock);

pc_form = 'nonwetting';
cap_scale = 10;
[kr, pc]  = tabulatedSatFunc([x, x.^2, y.^2, y.*cap_scale*barsa]);
%Define constant properties for viscosity and density

props = constantProperties([   1,  10] .* centi*poise, ...
                           [1000, 700] .* kilogram/meter^3);
%Here we put together a valid fluid object from the above defined functions. To read more about the fluid structure write help fluid_structure in MRST. First make a fluid without capillary pressure

fluid = struct('properties', props                  , ...
               'saturation', @(x, varargin)    x.s  , ...
               'relperm'   , kr);
%Then make another fluid object identical to the one above except for the capillary pressure term 'pc'.

fluid_pc = struct('properties', props                  , ...
                  'saturation', @(x, varargin)    x.s  , ...
                  'relperm'   , kr                     , ...
                  'pc'        , @(x, varargin) pc(x.s));
%Plot the pc-curve

%Make a dummy state/solution structure to plot the pc curve since 'fluid.pc' demands state as an input

xDummy   = initState(G, [], [0, 1]);
xDummy.s = linspace(0, 1, numel(xDummy.s))'; ...
pc = convertTo(fluid_pc.pc(xDummy), barsa);


figure
%clf
plot(xDummy.s, pc);
xlabel('s_w'); ylabel('pc [bar]');
title('Capillary pressure curve')
rate = 0.5*meter^3/day;
bhp  = 1*barsa;

%%
%Wells

W = verticalWell([], G, rock, 1, 1, 1:nz,          ...
                 'Type', 'rate', 'Val', rate, ...
                 'Radius', .1, 'Name', 'I', 'Comp_i', [1 0]);
W = verticalWell(W, G, rock, nx, ny, 1:nz,     ...
                 'Type','bhp', 'Val', bhp, ...
                 'Radius', .1, 'Dir', 'x', 'Name', 'P', 'Comp_i', [0 1]);
 figure(100)            
 plotWell(G, W);
 nw=length(W);
 %Z=eye(nx*nx);
 for nd=nx*ny+1:nx*ny+nw
    Z(nd,nd-nx*ny)=1;  
 end
 %Z=1;
 %figure
%spy(Z)
%Set up solution structures and assemble linear system

rSol    = initState(G, W, 0, [0.2, 0.8]);
rSol_pc = initState(G, W, 0, [0.2, 0.8]);

gravity off
verbose = false;

%S  = computeMimeticIP(G, rock, 'Verbose', verbose,'InnerProduct','ip_tpf');

%%
tol=1e-7;
  mrstModule add agmg
  %mrstModule add PCG_ICSol
   %  solver = AGMGSolverAD('tolerance', 1e-5);
   %solver = GMRES_ILUSolverAD('tolerance', 1e-5);
  
 %solver = PCG_ICSolverAD('tolerance', tol,'maxIterations', 1000);
 solver = DPCG_ICSolverAD('tolerance', tol,'maxIterations', 1000, 'Z',Z);
% solver = BackslashSolverAD();
% pressureSolver = BackslashSolverAD();
% linsolve = LinearSolverAD('ellipticSolver', pressureSolver);
%linsolver = LinearSolverAD('ellipticSolver', pressureSolver);

fn = @(A, b) solver.solveLinearSystem(A, b);



psolve  = @(state, fluid) incompTPFA(state, G, hT, fluid, 'wells', W,'MatrixOutput',true,'LinSolve', fn);

tsolve  = @(state, dT, fluid) implicitTransport(state, G, dT, rock, ...
                                                fluid, 'wells', W, ...
                                                'verbose', verbose);

rSol    = psolve(rSol, fluid);
rSol_pc = psolve(rSol_pc, fluid_pc);

T      = 300*day();
dT     = T/15;
dT = 300*day();
dTplot = 100*day();  % plot only every 100th day
N      = fix(T/dTplot);
pv     = poreVolume(G,rock);
%Start the main loop

t  = 0; plotNo = 1;
h1 = 'No pc - '; H2 = 'Linear pc - ';
e = []; p_org = []; p_pc = [];
figure;

while t < T,
    t
   % TRANSPORT SOLVE
   rSol    = tsolve(rSol, dT, fluid);
   rSol_pc = tsolve(rSol_pc, dT, fluid_pc);

   % Check for inconsistent saturations
   s = [rSol.s(:,1); rSol_pc.s(:,1)];
   assert(max(s) < 1+eps && min(s) > -eps);

   % Update solution of pressure equation.
   rSol    = psolve(rSol,    fluid);
   rSol_pc = psolve(rSol_pc, fluid_pc);

   % Measure water saturation in production cells in saturation
   e = [e; sum(abs(rSol.s(:,1) - rSol_pc.s(:,1)).*pv)/sum(pv)]; %#ok
   p_org = [p_org; rSol.s(W(2).cells,1)' ];                 %#ok
   p_pc  = [p_pc; rSol_pc.s(W(2).cells,1)'];                 %#ok

   % Increase time and continue if we do not want to plot saturations
   t = t + dT;
   if ( t < plotNo*dTplot && t<T) continue, end

   % Plot saturation
   heading = [num2str(convertTo(t,day)),  ' days'];
   r = 0.01;
   %subplot('position',[(plotNo-1)/N+r, 0.50, 1/N-2*r, 0.48]), cla
   figure
   plotCellData(G, rSol.s(:,1));
   caxis([0 1]), view(60,50), axis equal off, title([h1 heading])

   %subplot('position',[(plotNo-1)/N+r, 0.02, 1/N-2*r, 0.48]), cla
   figure
   plotCellData(G, rSol_pc.s(:,1));
   caxis([0 1]), view(60,50), axis equal off, title([H2 heading])

   plotNo = plotNo+1;

end
