% close all
% clear all
% rr=1;
% pod=0;
% cp=0;
% def=0;
%% Setting up the linear solver
%Select a linear solver 1.AGMG, 2.GMRES, 3.PCG_IC, 4.DPCG_IC, other Backslash
lsolver=3;
% Select transport solver 1. Explicit, 2. Implicit
tsolver=2;
% We define the maximum of iterations for the liner solver
maxIterations=500;
%steps =200;
% We define the tolerance of the linear solver
k=11;
tol = 10^-k;
% Compute estimated condition number
cn=0;
ceigs = 0;
        % Number of deflation vectors
dv=30;
podv{1}=[1:30];
podv{2}=[11:30];
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
gravity off

steps = 480;
cap_scale = 10;
muw = 0.5;
muo = 3;
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
%% Dimentions of the reservoir SPE 10
layers = 1;

nxi = 3;
nyi = 3;
nx = 16;
ny = 56;
nz = numel(layers);
Lx = 60;
Ly = 220;
Lz = 1;

    %The grid and the permeability are obtained from the MRST
   
    %physDims = [1200, 2200, 2*cartDims(end)] ;   % ft -> m

   % max(rock.perm)/min(rock.perm)
    %The permeability is 
    
 



%%
cartDims = [  nx,  ny, nz];
%cartDims = [nx ny nz];
domain = [Lx Ly Lz];
%physDims = [600, 2200, 2*cartDim(end)] ;  
physDims = [1200, 4200, 2*cartDims(end)] .* ft();
%physDims = cartDims .* [20, 10, 2]*ft;
%G = cartGrid([nx ny nz]);
%G = computeGeometry(G);
G      = computeGeometry(cartGrid(cartDims,physDims));
%G = computeGeometry(cartGrid(cartDims, physDims));

rock     = SPE10_rock(layers);
rock.perm = convertFrom(rock.perm, milli*darcy);
    perm = rock.perm(:,1);
    perm = reshape(perm,domain);
    permupscaled = sampleFromBox(G,perm);
    rock.perm=permupscaled;
   poro = rock.poro(:,1);
    poro = reshape(poro,domain);
    poroupscaled = sampleFromBox(G,poro);
    rock.poro=poroupscaled;    
    
    
%     max(rock.perm)
%     min(rock.perm)
    contperm=max(rock.perm)./min(rock.perm);
is_pos             = rock.poro > 0;
rock.poro(~is_pos) = min(rock.poro(is_pos));



%G      = computeGeometry(cartGrid(cartDim,physDims));

%% Wells
%W1_val = 5*meter^3/day;
% W1_val = 100*barsa;
% W2_val = 0*barsa;
% %% Wells values
well(2:5)=210*barsa;
% % well(5)=5*meter^3/day;
well(1) = 700*barsa;
well(6) = 700*barsa;
%% Rock properties

% Set rock initial properties
%rock   = makeRock(G, rperm*milli*darcy, rporo);
pv     = poreVolume(G, rock);




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
%wrad     = [0.125, 0.125, 0.125, 0.125, 0.125] .* meter;
wrad(1:6) = 1*meter;
%{Sub2ind(ceil(nx/2),ceil(ny/2),1:nz,nx,ny,nz),......
wloc     =  {Sub2ind(ceil(nxi),ceil(ny/2-5),1:nz,nx,ny,nz), ...
            Sub2ind(nxi,nyi,1:nz,nx,ny,nz),Sub2ind(nxi,ny-5,1:nz,nx,ny,nz), ...
            Sub2ind(nx-3,nyi,1:nz,nx,ny,nz),Sub2ind(nx-3,ny-5,1:nz,nx,ny,nz)};
wname    = {'I1', 'P1', 'P2', 'P3', 'P4'};
Sign      = [ 1 ,  -1 ,  -1 ,  -1 , -1 ];
Comp_iI= [1 0];
Comp_iP= [0 1];
% 6th well
wtype{6} = 'bhp';
wtarget(6) = well(6);
%wrad(6) = 0.125*meter;
wloc{6} = Sub2ind(ceil(nx-3),ceil(ny/2+5),1:nz,nx,ny,nz);
wname{6} = 'I2';
Sign(6) = 1;


for i = 1 : numel(wtype)
if Sign(i) == -1
 wcomp{i} = Comp_iP;
else
wcomp{i} = Comp_iI;
end
end

W = [];
for w = [1 2 3 4 5 6]
   W = addWell(W, G, rock, wloc{w}, ...
                    'Type', wtype{w}, 'Val', wtarget(w), ...
                    'Radius', wrad(w), 'Name', wname{w}, ...
                    'Comp_i', wcomp{w}, ...
                    'Sign', Sign(w));
end

% 
% nf = nf + 1;
% f(nf) = figure(nf);
% figure(nf)
% %plotCellData(G, log(rock.perm(:,1)),'LineStyle','none');
% 
% clf;
% plotCellData(G, log(rock.perm(:,1)));
% for i = 2 :numel(W)
%     plotWell(G, W(i),'color', 'r');
% end
% plotWell(G, W(1), 'color', 'b');
% file{nf} = 'Permeability';
% %title([file{nf} ' field [mD]'])
% if nz > 1
%     view(50,30)
% else
%     view(0,90)
% end
% axis equal off
% hold off
% break




close all
nf = nf + 1;
f(nf) = figure(nf);
figure(nf)
h=plotCellData(G, log(rock.perm(:,1)),'LineStyle','none'); 
axis equal tight off
% xmi = min(log(rock.perm(:,1)));
% xma = max(log(rock.perm(:,1)));
% clim = [xmi xma ];
% 
% caxis(clim);
% %
% h=colorbar;
% set(h,'CLim',clim);
%  set(h,'position',[0.65 0.2 0.05 0.70])
for i = 2 :5
    plotWell(G, W(i),'color', 'r');
end
plotWell(G, W(1), 'color', 'b');
plotWell(G, W(6), 'color', 'b');
file{nf} = 'Permeability';
%title([file{nf} ' field [mD]'],'position',[0.45 0.05 0.35 0.40])
if nz > 1
    view(10,20)
else
    view(0,90)
end
hold off
%break

                    %% Changing wells
                    %Create a well structure to support multiple report steps
                    %clear newW
                    [newW{1:steps,1}] = deal(W);
                     BHP = zeros(numel(W),steps);
                     
%                    % Valr = zeros(numel(W),steps);
%                     %Then we'll assign the time-dependent controls
                    tch =5;
                    csteps = round(steps/tch);
%                     
                    h = 0;
                    while tch*(h+1) < steps
                        h=h+1;
                       BHP(1,(h-1)*tch+1:tch*h) = rand;
                       BHP(3,(h-1)*tch+1:tch*h) = rand;
                    end
% %                     
                     BHP(2,:) = 1-BHP(3,:);
                     BHP(4,:) = 1-BHP(1,:);
%                     figure
%                     plot(time,BHP(1,:))
%                     hold on
%                     plot(time,BHP(2,:))
%                     hold on
%                     plot(time,BHP(3,:))
%                     hold on
%                     plot(time,BHP(4,:))
%                     axis tight
%                     break
%                    
                   
%                     Valr(1,1:steps/2) = linspace(0.5, 0.8, steps/2)/day;
%                     Valr(2,steps/2+1:steps) = linspace(0.5, 0.8, steps/2)/day;
%                     Valr(4,1:steps/2) = linspace(0.5, 0.8, steps/2)/day;
%                     Valr(3,steps/2+1:steps) = linspace(0.5, 0.8, steps/2)/day;
                    %Valr(5,1:steps) = linspace(0.5, 0.8, steps)/day;
                    %Valr =  linspace(0.5, 1, steps)/day();
                    W1=newW;
                    for w = 1:numel(newW),
                        for i = 2:5
                            newW{w}(i).type = 'bhp';
                             newW{w}(i).val= (BHP(i-1,w)*140+140)*barsa;
%                             newW{w}(i).type = 'rate';
%                            newW{w}(i).val = Valr(i,w);                       
                        end
                    end
                    
                    W0=newW;







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
h=plot(s, kr); legend('kr_1', 'kr_2'),title(['Relative permeability']), axis equal tight
hold off

%% Create the directory

%dir='/mnt/sda2/cortes/Results/17_08/31/1/';
dir='/mnt/sda2/cortes/Results/2017/Report/wbt/SPE10/training/1/';
folder=[ num2str(nx) '_' num2str(ny) '_l_' num2str(numel(layers)) 'cp_'  num2str(cp) ];
mkdir([dir], folder)
dir1 = [dir folder '/'];

folder=[  'def_' num2str(def) '_pod_' num2str(numel(dpod)) ];
mkdir([dir1], folder)
dir2 = [dir1 folder '/'];
if training == 0

folder=[  't_' num2str(training) ];
mkdir([dir1], folder)
dir2 = [dir1 folder '/'];
end



pmark =['o' '+' '*' 'x' 'd' 's' ];
pcol ={'r' [0 0.7 0] 'b'  [0 0.7 0.7] [0.7 0 0.7] [0.5 0.5 0.7] };