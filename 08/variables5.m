%% Setting up the linear solver
%Select a linear solver 1.AGMG, 2.GMRES, 3.PCG_IC, 4.DPCG_IC, other Backslash
lsolver=3;
% Select transport solver 1. Explicit, 2. Implicit
tsolver=2;
% We define the maximum of iterations for the liner solver
maxIterations=500;
% We define the tolerance of the linear solver
k=11;
tol = 10^-k;
% Compute estimated condition number
cn=0;
ceigs = 0;
        % Number of deflation vectors
dv=10;
podv{1}=[1:10];
podv{2}=[6:10];
dpod=podv{rr};
    if pod== 0
                dpod = [];
            end 
            
%% Save figures option
sf = true;
% Numbering for figures
nf =0;
%% Set up model
% Square domain with homogeneous rock properties, no flow across the
% boundaries, no gravity, injection in the southeast and production in the
% northwest corners, both operating at fixed bottom-hole pressure
gravity reset off


cap_scale = 10;
muw = 1;
muo = 10;
rhow = 1000;
rhoo = 700;
Sro = 0;
Srw = 0;
x = linspace(0, 1, 11) .';
y = linspace(1, 0, 11) .';
nw = 2;
no = 2;
kw = 1;
ko = 1;
%% Dimentions of the reservoir
sz = 35;
nxi=1;
nyi=1;
nx=sz;
ny=sz;
nz=1;
Lx=sz;
Ly=sz;
Lz = 1;

%%
cartDim = [nx ny nz];
domain = [Lx Ly Lz];
G      = computeGeometry(cartGrid(cartDim,domain));

%% Rock properties
rperm = 100;
rporo = 0.25;
%  Number of layers with different permeability
pl = true;
rlay = 5;
%per = 0;

%% Wells
%W1_val = 5*meter^3/day;
% W1_val = 100*barsa;
% W2_val = 0*barsa;
% %% Wells values
well(2:5)=-50*barsa;
% % well(5)=5*meter^3/day;
well(1) = 200*barsa;


%% Rock properties

% Set rock initial properties
rock   = makeRock(G, rperm*milli*darcy, rporo);
pv     = poreVolume(G, rock);

% Create layers of  diverse permeability
if pl == true
    v = [];
    for i = 1:2*rlay:ny
        for j = 0:rlay-1
            if i+j < ny+1
                v = [v j+i];
            end
        end
    end


    [I] = Sub2ind([1:nx],v,1:nz,nx,ny,nz);
    rock.perm(I) = rperm*10^(-per)*milli*darcy();

end




%% Add 2 wells
% [I1] = Sub2ind(1,1,1:nz,nx,ny,nz);
% [I2] = Sub2ind(nx,ny,1:nz,nx,ny,nz);
% W = addWell([],G, rock, I1, 'Type', 'bhp', ...
%     'Val', W1_val , 'name', 'I', 'radius', .1, 'Comp_i', [1 0]);
% W = addWell(W,G, rock, I2, 'Type', 'bhp', ...
%     'Val', W2_val , 'name', 'P', 'radius', .1, 'Comp_i', [0 1]);

%% Various wells
wtype    = {'bhp', 'bhp', 'bhp', 'bhp', 'bhp'};
wtarget  = [well(1),   well(2),   well(3),   well(4), well(5)];
wrad     = [0.125, 0.125, 0.125, 0.125, 0.125] .* meter;
wloc     = {Sub2ind(ceil(nx/2),ceil(ny/2),1:nz,nx,ny,nz),......
            Sub2ind(nxi,nyi,1:nz,nx,ny,nz),Sub2ind(nxi,ny,1:nz,nx,ny,nz), ...
            Sub2ind(nx,nyi,1:nz,nx,ny,nz),Sub2ind(nx,ny,1:nz,nx,ny,nz)};
wname    = {'I', 'P1', 'P2', 'P3', 'P4'};
Sign      = [ 1 ,  -1 ,  -1 ,  -1 , -1 ];
Comp_iI= [1 0];
Comp_iP= [0 1];

for i = 1 : numel(wtype)
if Sign(i) == -1
 wcomp{i} = Comp_iP;
else
wcomp{i} = Comp_iI;
end
end

W = [];
for w = [1 2 3 4 5]
   W = addWell(W, G, rock, wloc{w}, ...
                    'Type', wtype{w}, 'Val', wtarget(w), ...
                    'Radius', wrad(w), 'Name', wname{w}, ...
                    'Comp_i', wcomp{w}, ...
                    'Sign', Sign(w));
end




nf = nf + 1;
f(nf) = figure(nf);
figure(nf)
plotCellData(G, rock.perm/(milli*darcy()),'LineStyle','none'); colorbar
 axis equal tight off
for i = 2 :numel(W)
    plotWell(G, W(i),'color', 'r');
end
plotWell(G, W(1), 'color', 'b');
file{nf} = 'Permeability';
title([file{nf} ' field [mD]'])
if nz > 1
    view(10,20)
else
    view(0,90)
end
hold off

%% Fluid model

%% Fluid with no capillary pressure
% fluid = initCoreyFluid('mu' , [ muw,  muo]*centi*poise     , ...
%     'rho', [rhow, rhoo]*kilogram/meter^3, ...
%     'n'  , [   nw,  no]                 , ...
%     'sr' , [  Srw, Sro]                 , ...
%     'kwm', [   kw, ko]);
% s=linspace(0,1,50)';
% nf = nf + 1;
% f(nf) = figure(nf);
% file{nf} = ['Relative_Permeability'];
% clf
% [kr,dkr,ddkr] = fluid.relperm(s);
% plot(s(s<=1-Sro),kr(s<=1-Sro,1),s(s>=Srw),kr(s>=Srw,2),'LineWidth',1); axis([0 1 -.02 1]);
% legend('kr_w', 'kr_o'), axis equal tight

% Set up capillarity and relative permeability

pc_form = 'nonwetting';
[kr, pc]  = tabulatedSatFunc([x, x.^2, y.^2, y.*cap_scale*barsa]);
props = constantProperties([   muw,  muo] .* centi*poise, ...
    [rhow, rhoo] .* kilogram/meter^3);


if cp==0
    fluid = initSimpleFluid('mu' , [   muw,    muo] .* centi*poise     , ...
        'rho', [rhow, rhoo] .* kilogram/meter^3, ...
        'n'  , [   nw,    no]);
else
    
    fluid = struct('properties', props                  , ...
        'saturation', @(x, varargin)    x.s  , ...
        'relperm'   , kr                     , ...
        'pc'        , @(x, varargin) pc(x.s));
    xDummy   = initState(G, [], [0, 1]);
    xDummy.s = linspace(0, 1, numel(xDummy.s))'; ...
        pc = convertTo(fluid.pc(xDummy), barsa);
    
    nf = nf + 1;
    figure(nf)
    file{nf} = ['Capillary_Pressure'];
    figure(nf)
    h=plot(xDummy.s, pc);
    xlabel('S_w'); ylabel('Pc [bar]');
    title('Capillary pressure curve')
    hold off
end
s=linspace(0,1,50)';


nf = nf + 1;
figure(nf)
file{nf} = ['Relative_permeability'];
s = linspace(0, 1, 1001).'; kr = fluid.relperm(s);
h=plot(s, kr); %legend('kr_1', 'kr_2')%,title(['Relative permeability']), axis equal tight
hold off

%% Create the directory

dir='/mnt/sda2/cortes/Results/2017/Report/5wells1/';

folder=[ '10-' num2str(k) '_' num2str(sz) 'perm_' num2str(per) 'cp' num2str(cp)];
mkdir([dir], folder)
dir1 = [dir folder '/'];

folder=[  'def_' num2str(def) '_pod_' num2str(numel(dpod)) ];
mkdir([dir1], folder)
dir2 = [dir1 folder '/'];


pmark =['o' '+' '*' 'x' 'd' ];
pcol ={'r' [0 0.7 0] 'b'  [0 0.7 0.7] [0.7 0 0.7] };
