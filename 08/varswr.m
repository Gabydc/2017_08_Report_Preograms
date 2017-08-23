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
nz=10;
% Compute estimated condition number
cn=0;
ceigs = 0;
%layers = 1:85;

Lx=300;
Ly=300;
Lz = 10;
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
muo = 1;
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
% Number of deflation vectors
dv=10;
podv{1}=[1:10];
podv{2}=[6:10];
%% Transport loop vars
T      = 100*day();
%T      = 100*day();
steps = 1000;
dT     = T/steps;
t = 0 : dT : T;
dTplot = ceil(T/3);  % plot only every 100th day
