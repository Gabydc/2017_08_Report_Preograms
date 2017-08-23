
cp = [0 ];

% The contrast between permeability layers can be varied one layer has
% permeability of 10*milli*darcy, the secon is varied as
% 10*10^(-per)*milli*darcy
per=[0];

% In this part the solver is chosen 0 means ICCG (no deflation) and
% 1 is with deflation
def=[0 ];

% If we want to use POD as deflation vectors pod=1
pod = [0];
dpod =0;
% Initialize the variables, size, fluid values, time
% step, ...
varswr
close all

%Create the directory
dir='/mnt/sda2/cortes/Results/17_05/two_phases/eigs/28/';

folder=[ '10-' num2str(k) '_' num2str(sz) 'nz' num2str(nz) 'perm_' num2str(per) 'cp' num2str(cp)];
mkdir([dir], folder)
dir1 = [dir folder '/'];

folder=[  'def_' num2str(def) '_pod_' num2str(numel(dpod)) ];
mkdir([dir1], folder)
dir2 = [dir1 folder '/'];

%Create the grid
G = cartGrid([sz, sz, nz], [sz, sz, nz]);
G = computeGeometry(G);

% Set up uniform permeability and constant porosity
rock.perm = ones(G.cells.num, 1)*10*milli*darcy();
rock.poro = ones(G.cells.num, 1)*0.2;

%                     % Create layers of permeability
%                     nlay = 1;
%                     v = [];
%                     for i = 1:2*nlay:ny
%                         for j = 0:nlay-1
%                             if i+j < ny+1
%                                 v = [v j+i];
%                             end
%                         end
%                     end
%
%
%                     [I] = Sub2ind([1:nx],v,1:nz,nx,ny,nz);
%                     rock.perm(I) = 10*10^(-per)*milli*darcy();


%% Compute half transmissibilities
% All we need to know to develop the spatial discretization is the reservoir
% geometry and the petrophysical properties. This means that we can compute
% the half transmissibilities without knowing any details about the fluid
% properties and the boundary conditions and/or sources/sinks that will
% drive the global flow:
hT = simpleComputeTrans(G, rock);

%% Fluid model
% We set up a two-phase fluid. Viscosity, density is
% set for oil and water using the valuers from vars

verbose = false;
if cp==0
    fluid = initSimpleFluid('mu' , [   muw,    muo] .* centi*poise     , ...
        'rho', [rhow, rhoo] .* kilogram/meter^3, ...
        'n'  , [   krw,    kro]);
else
    
    fluid = struct('properties', props                  , ...
        'saturation', @(x, varargin)    x.s  , ...
        'relperm'   , kr                     , ...
        'pc'        , @(x, varargin) pc(x.s));
    xDummy   = initState(G, [], [0, 1]);
    xDummy.s = linspace(0, 1, numel(xDummy.s))'; ...
        pc = convertTo(fluid.pc(xDummy), barsa);
    
    hcp=figure(232);
    plot(xDummy.s, pc);
    xlabel('s_w'); ylabel('pc [bar]');
    title('Capillary pressure curve')
    
end

% Plot relative permeability
s = linspace(0, 1, 1001).'; kr = fluid.relperm(s);
f(2) = figure(2);
plot(s, kr); legend('kr_1', 'kr_2'),title(['rel perm']), axis equal tight
%Boundary conditions. We set boundary conditions,
%waterflooding, right boundary / wells
pv = poreVolume(G, rock);
%injRate = -sum(pv)/(500*day);
%                     injRate = -0.1*meter^3/day;
%                     bc = fluxside([], G, 'xmin', -injRate, 'sat', [1, 0]);
%                     bc = pside(bc, G, 'xmax', 0*barsa, 'sat', [0, 1]);
%
W = addWell([],G, rock, 1:nx*ny:nx*ny*(nz), 'Type', 'rate', ...
    'Val', 0.5/day(), 'name', 'I', 'radius', .1, 'Comp_i', [1 0]);
W = addWell(W,G, rock, nx*ny:nx*ny:G.cells.num, 'Type', 'bhp', ...
    'Val', 0, 'name', 'P', 'radius', .1, 'Comp_i', [0 1]);
f(1)=figure(1);
plotCellData(G, rock.perm/(milli*darcy())); colorbar
colormap(jet), axis equal tight off
plotWell(G, W(1),'color', 'r');
plotWell(G, W(2), 'color', 'b');
title('Permeability field [mD]')
if nz > 1
    view(10,20)
else
    view(0,90)
end

% Create a well structure to support multiple report steps
[newW{1:steps}] = deal(W);

%Then we'll assign the time-dependent controls
BHP = 100*rand(steps,2);
for w = 1:numel(newW),
for i = 1:numel(W)
    newW{w}(i).type = 'bhp';
    newW{w}(i).val = BHP(w,i);
  end  
end

% Adjust the constraints structure to account for the time-dependent controls

% 
% [schedule.control(1:steps).W] = newW{:};
% 
% schedule.step.control = 1:numel(newW);
% 
% schedule.step.val = diff([zeros(size(t),1) ; t]);
% break