%% Set up grid and petrophysical data
% We use a Cartesian grid of size nx-by-ny with homogeneous petrophysical
% data: permeability of 100 mD and porosity of 0.2.
dir='/mnt/sda2/cortes/Results/17_06/two_phases/eigs/17/';
%Select a linear solver 1.AGMG, 2.GMRES, 3.PCG_IC, 4.DPCG_IC, other Backslash
lsolver=3;
% Select transport solver 1. Explicit, 2. Implicit
tsolver=2;
% We define the size of the reservoir and the inicial points of the grid
% sz is the number of cells, is the same for the two dimensions
sz=30;


% If we use wells, we can define here the position of the wells
% w1=10*2^pw;
% w2=2*w1;
nxi=1;
nyi=1;
nx=sz;
ny=sz;
nz=1;
% Compute estimated condition number
cn=0;
ceigs = 0;
%layers = 1:85;

Lx=300;
Ly=300;
Lz = 1;
gravity off

%% Create the grid
G = cartGrid([sz, sz, nz], [Lx, Ly, Lz]);
G = computeGeometry(G);

%% Set up uniform permeability and constant porosity
per1 = 100;
por1 = 0.25;
rock.perm = ones(G.cells.num, 1)*per1*milli*darcy();
rock.poro = ones(G.cells.num, 1)*por1;

% We define the maximum of iterations for the liner solver
maxIterations=500;
% We define the tolerance of the linear solver
k=7;
tol = 10^-k;
% We can change the permeability one layer is 10*milli*darcy(), the second
% is 10*milli*darcy()*10^(-per)


pc_form = 'nonwetting';
muw = 1;
muo = 10;
rhow = 1000;
rhoo = 700;
Smaxo = 0.5;
Smaxw = 0.5;
Sro = 0.2;
Srw = 0.2;
x = linspace(0, 1, 11) .';
y = linspace(1, 0, 11) .';
krw = 1.5;
kro = 2;
krw0 = 0.4;
kro0 = 1;
cap_scale = 10;
props = constantProperties([   muw,  muo] .* centi*poise, ...
    [rhow, rhoo] .* kilogram/meter^3);
[kr, pc]  = tabulatedSatFunc([x, x.^krw, y.^kro, y.*cap_scale*barsa]);
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
                      well(1:4)=76 * barsa;
                   well(5)=5*meter^3/day;  
                  % well(5) = 0.3;
                  
wtype    = {'bhp', 'bhp', 'bhp', 'bhp', 'rate'};
wtarget  = [well(1),   well(2),   well(3),   well(4), well(5)] ;
wrad     = [0.125, 0.125, 0.125, 0.125, 0.125] .* meter;
wloc     = [  nxi,   nxi,     nx,   nx, nx;
              nyi,   ny,     nyi,   ny, ny];
%wname    = {'I1', 'I2', 'I3', 'I4', 'P'};
wname    = {'P1', 'P2', 'P3', 'P4', 'I'};
sgn      = [ -1 ,  -1 ,  -1 ,  -1 , 1 ];
Comp_iI= [1 0];
Comp_iP= [0 1];
%%
for i = 1 : numel(wtype)
if sign(i) == -1
 wcomp{i} = Comp_iP;
else
wcomp{i} = Comp_iI;
end
end

W = [];        
for w = 1 : numel(wtype),
   W = verticalWell(W, G, rock, wloc(1,w), wloc(2,w), [], ...
                    'Type', wtype{w}, 'Val', wtarget(w), ...
                    'Radius', wrad(w), 'Name', wname{w}, ...
                    'Comp_i', wcomp{w}, ...
                    'Sign', sgn(w), 'InnerProduct', 'ip_tpf');
end

%                     W = addWell([],G, rock, 1:nx*ny:nx*ny*(nz), 'Type', 'rate', ...
%                         'Val', 0.5/day(), 'name', 'I', 'radius', .1, 'Comp_i', [1 0]);
%                     W = addWell(W,G, rock, nx*ny:nx*ny:G.cells.num, 'Type', 'bhp', ...
%                         'Val', 0, 'name', 'P', 'radius', .1, 'Comp_i', [0 1]);
  
                    f(1)=figure(1);
                    plotCellData(G, rock.perm/(milli*darcy())); colorbar
                    colormap(jet), axis equal tight off
for i = 1 :numel(W)-1
                    plotWell(G, W(i),'color', 'r');
end
                    plotWell(G, W(numel(W)), 'color', 'b');
                    title('Permeability field [mD]')
                    if nz > 1
                        view(10,20)
                    else
                        view(0,90)
                    end
                    
% Number of deflation vectors
dv=10;
podv{1}=[1:10];
podv{2}=[6:10];


%% Initial pressure
pinit = 0; 
%% Transport loop vars
T      = 1000*day();
%T      = 100*day();
steps = 100;
dT     = T/steps;
t = 0 : dT : T;
dTplot = ceil(T/3);  % plot only every 100th day

