%% Simulate the SPE10 base case
% The simulation uses a mimetic pressure solver and an implicit transport
% solver.

clear, close all hidden

%%
use_mimetic = false;
use_reorder = false;
use_ICCG = true;
use_DICCG = false;
spe10_data  = fullfile(fileparts(mfilename('fullpath')), ...
                       '..', 'spe10_rock.mat');
if ~exist(spe10_data, 'file'),
   if ~make_spe10_data,
      error(['Failed to establish on-disk representation of ', ...
             'SPE10 rock data']);
   end
end

%%
layers = 1 : 35;

cartDims = [60, 220, numel(layers)];
rock     = getSPE10rock(layers);

rock.perm = convertFrom(rock.perm, milli*darcy);

is_pos             = rock.poro > 0;
rock.poro(~is_pos) = min(rock.poro(is_pos));

physDims = cartDims .* [20, 10, 2]*ft;

G = computeGeometry(cartGrid(cartDims, physDims));

if use_mimetic,
   try
      require mimetic
   catch
      mrstModule add mimetic
   end

   S = computeMimeticIP(G, rock);
else
   T = computeTrans(G, rock);
end

%%
fluid = initSimpleFluid('mu' , [   1,  10]*centi*poise     , ...
                        'rho', [1014, 859]*kilogram/meter^3, ...
                        'n'  , [   2,   2]);

%%
% Set Comp_i=[0,0] in producers to counter X-flow effects...
%
% Initial pressure in the producers
ipp = 400;
% Initial rate in the injector
iir = 5000;
if use_mimetic,
   well_ip = 'ip_simple';
else
   well_ip = 'ip_tpf';
end
W = verticalWell([], G, rock,  1,   1, [], 'Type', 'bhp', ...
                 'InnerProduct', well_ip, ...
                 'Val', ipp*psia, 'Radius', 0.125*meter, ...
                 'Name', 'P1', 'Comp_i', [0, 0]);

W = verticalWell(W , G, rock, 60,   1, [], 'Type', 'bhp', ...
                 'InnerProduct', well_ip, ...
                 'Val', ipp*psia, 'Radius', 0.125*meter, ...
                 'Name', 'P2', 'Comp_i', [0, 0]);

W = verticalWell(W , G, rock, 60, 220, [], 'Type', 'bhp', ...
                 'InnerProduct', well_ip, ...
                 'Val', ipp*psia, 'Radius', 0.125*meter, ...
                 'Name', 'P3', 'Comp_i', [0, 0]);

W = verticalWell(W , G, rock,  1, 220, [], 'Type', 'bhp', ...
                 'InnerProduct', well_ip, ...
                 'Val', ipp*psia, 'Radius', 0.125*meter, ...
                 'Name', 'P4', 'Comp_i', [0, 0]);

W = verticalWell(W , G, rock, 30, 110, [], 'Type', 'rate',   ...
                 'InnerProduct', well_ip, ...
                 'Val', iir*stb/day, 'Radius', 0.125*meter, ...
                 'Name', 'I1', 'Comp_i', [1, 0]);

%% Set time
DT    = 100*day;
nstep =  5;      % 2000 days

%% Changing wells
%Create a well structure to support multiple report steps
%clear newW
[newW{1:nstep,1}] = deal(W);
BHP = zeros(numel(W),nstep);

%                    % Valr = zeros(numel(W),nstep);
%                     %Then we'll assign the time-dependent controls
tch =10;
csteps = round(nstep/tch);
%
h = 0;
while tch*(h+1) < nstep
    h=h+1;
    BHP(1,(h-1)*tch+1:tch*h) = rand;
    BHP(3,(h-1)*tch+1:tch*h) = rand;
end
%
BHP(2,:) = 1-BHP(3,:);
BHP(4,:) = 1-BHP(1,:);
W1=newW;
for w = 1:numel(newW),
    for i = 1:4
        newW{w}(i).type = 'bhp';
        newW{w}(i).val= (BHP(i,w)*ipp)*psia;
    end
end

W0=newW;
% 
% time = 0 : DT : DT*(nstep-1);
% 
%                     figure
%                     for i = 1 : numel(W)-1
%                         for j = 1 : nstep
%                         wval(j)=W0{j}(i).val;
%                         end
%                     plot(time,wval)
%                     hold on
%                     end
%                     figure
%                     plot(time,BHP(1,:))
%                     hold on
%                     plot(time,BHP(2,:))
%                     hold on
%                     plot(time,BHP(3,:))
%                     hold on
%                     plot(time,BHP(4,:))
%                     axis tight
%                    % break
% 


             
             
%%
x         = initResSol (G, 0);
x.wellSol = initWellSol(W, 0);

%%
tol = 5.0e-11;
maxIterations = 1000;
cn = 0;

p0 = x.pressure;
n = size(p0);
for i=1:numel(W)
    p0(n+i) = 0;
end


linsolve_t = @(J, F) agmg(J, F, 50, 5.0e-11, 2000, 0);


if ~use_reorder,
   tsolve = @(x, dt) ...
      implicitTransport(x, G, dt, rock, fluid, 'wells', W, ...
                        'LinSolve', linsolve_t);
else
   mu          = fluid.properties();
   fluid.param = struct('viscw', mu(1), 'visco', mu(2),       ...
                        'srw', 0, 'sro', 0, 'nw', 1, 'no', 1, ...
                        'satnum', ones([G.cells.num, 1], 'int32'));
   tsolve = @(x, dt) ...
      implicitTransportReorder(x, G, dt, rock, fluid, 'wells', W);
end

%%


Prod = struct('t'  , []                  , ...
              'vpt', zeros([0, numel(W)]), ...
              'opr', zeros([0, numel(W)]), ...
              'wpr', zeros([0, numel(W)]), ...
              'wc' , zeros([0, numel(W)]));

append_wres = @(x, t, vpt, opr, wpr, wc) ...
   struct('t'  , [x.t  ; t                  ], ...
          'vpt', [x.vpt; reshape(vpt, 1, [])], ...
          'opr', [x.opr; reshape(opr, 1, [])], ...
          'wpr', [x.wpr; reshape(wpr, 1, [])], ...
          'wc' , [x.wc ; reshape(wc , 1, [])]);

wres = cell([1, 4]);

for k = 1 : nstep,
    W = W0{k};
 p0(1:n) = x.pressure;
 if use_ICCG,
solver = PCG_ICSolverAD_cn('tolerance', tol,'maxIterations', maxIterations,'cn',cn,'x0',p0);
linsolve_p = @(A, b) solver.solveLinearSystem(A, b);
    else
linsolve_p = @(S, h) agmg(S, h,  1,  tol,  maxIterations, 0);
end
if use_mimetic,
   psolve = @(x) ...
      incompMimetic(x, G, S, fluid, 'wells', W, 'LinSolve', linsolve_p);
else if use_ICCG,
  
   psolve = @(x,p0) ...
      incompTPFA_g_o(x, G, T, fluid, 'wells', W,'LinSolve', linsolve_p);
    else
   psolve = @(x) ...
      incompTPFA_amg(x, G, T, fluid, 'wells', W, 'LinSolve', linsolve_p);
    end
end

 %   if use_ICCG,
   t0 = tic; [x, report(k)]  = psolve(x,p0);     dt = toc(t0);
  %  else
    % t0 = tic;  x =  psolve(x);     dt = toc(t0);
  %  end
   fprintf('[%02d]: Pressure:  %12.5f [s]\n', k, dt);

   t0 = tic; x = tsolve(x, DT); dt = toc(t0);
   fprintf('[%02d]: Transport: %12.5f [s]\n', k, dt);

   t = k * DT;

   [wres{:}] = prodCurves(W, x, fluid);
   Prod(k)      = append_wres(Prod, t, wres{:});
   Pressure(k) = x.pressure;
end


%%

                if training == 0
                    
                    Z =  Pressure;
                    %[U,S]=defpodf_Dt(Z,dir2,dv,ts,dTplot/day());
                    [U,S]=PODbasis(Z);
                    
                    nf = nf + 1;
                    file{nf} = ['eig_pod_' num2str(t/day())];
                    f(nf) = figure(nf);
                    plot(log((diag(S))),'*r');
                    ylabel('log(Value) ','FontSize',16)
                    xlabel('Eigenvalue','FontSize',16)
                    axis('tight')
                    Z=U(:,dpod);  
                    for i=1:numel(W)
                    
                    Z(nx*ny+i,:)=0;
                    end
                    save([dir1 'Z.mat'],'Z')
                    
                    save([dir1 'T.mat'],'T')
                    
                end








%%
figure
plotCellData(G, x.s(:,1), 'EdgeColor', 'k', ...
             'EdgeAlpha', 0.050, 'FaceAlpha', 0.375)
view(3), colorbar, axis tight off

figure
plotCellData(G, x.s(:,1), find(x.s > 0.5), 'EdgeColor', 'k', ...
             'EdgeAlpha', 0.050, 'FaceAlpha', 0.375)
view(3), colorbar, axis tight off

figure
plot(convertTo(Prod.t, day), convertTo(Prod.vpt(:,1:end-1), meter^3/day))
legend({ W(1:end-1).name }, 'Location', 'Best')
xlabel('Time [d]'), ylabel('Total Production Rate [m^3/d]')

figure
plot(convertTo(Prod.t, day), convertTo(Prod.opr(:,1:end-1), meter^3/day))
legend({ W(1:end-1).name }, 'Location', 'Best')
xlabel('Time [d]'), ylabel('Oil Production Rate [m^3/d]')

figure
plot(convertTo(Prod.t, day), convertTo(Prod.wpr(:,1:end-1), meter^3/day))
legend({ W(1:end-1).name }, 'Location', 'Best')
xlabel('Time [d]'), ylabel('Water Production Rate [m^3/d]')

figure
plot(convertTo(Prod.t,day), Prod.wc(:,1:end-1))
legend({ W(1:end-1).name }, 'Location', 'Best')
xlabel('Time [d]'), ylabel('Well Water Cut')




% 
%                 Np = 3;
%                 time = (0:DT:t)/day;
% %                 px = [0.01 0.26 0.51 0.76];
% %                 nz=35;
% %                 nf = 11;
%                  figure
% 
%                  plot(report(1).residual)
% %                 %title('Water Saturation');
% %                 file{nf} = ['Water_saturation'];
% %                 clim = [0 1];
% %                 subplotcbspe(nf,clim,t,Np,G,nz,T,x.s)
% % 
% 
% 




