
    close all
clear all
clc
varsspe
    cartDims = [  nx,  ny,   numel(layers)];
    %physDims = [nx, ny, 2*cartDims(end)] ;   % ft -> m
    G = cartGrid([nx ny numel(layers)]);
    G = computeGeometry(G);

    % Set up uniform permeability and constant porosity
       rock = getSPE10rock(layers);
low = 1e-4;
rock.poro(rock.poro < low) = low;
    %The permeability is 
    perm = rock.perm(:,1);
 


%%
for cp = [0 1]
   

   
for def=[0  1]
  
    
    

% Number of deflation vectors

% If we want to use POD pod==1
for pod = [0 1] 
    clearvars -except per def pod cp cartDims G rock perm


dpod=[];
if (def == 0)  && (pod == 1)
           pod=0;
    continue
end
     for rr=1:2
         if (pod == 0) && (rr==2) 
             continue 
         end
         varsspe
         close all
                dpod=podv{rr};
                if pod== 0
                    dpod=[];
                end
                
     f(1)=figure(1);
      
            plotCellData(G, log(perm),'linestyle','none');
            view(50,15) 
            colormap(jet) 
            axis equal tight 
            set(gca,'position',[0.01 0.01 0.9 0.9])
            
    
    %Create the directory
    dir='/mnt/sda2/cortes/Results/17_08/spe/';
    
    folder=[ 'spe10-' num2str(k) '_' 'nz' num2str(nz) 'cp' num2str(cp)];
    mkdir([dir], folder)
    dir1 = [dir folder '/'];
        
    folder=[  'def_' num2str(def) '_pod_' num2str(numel(dpod)) ];
    mkdir([dir1], folder)
    dir2 = [dir1 folder '/'];
    %handler2 = ResultHandler('writeToDisk', true,'dataDirectory',dir,'dataFolder','results');
    %Create the grid
    %domain = [250 250 20];


   % break
    %% Compute half transmissibilities
    % All we need to know to develop the spatial discretization is the reservoir
    % geometry and the petrophysical properties. This means that we can compute
    % the half transmissibilities without knowing any details about the fluid
    % properties and the boundary conditions and/or sources/sinks that will
    % drive the global flow:
    hT = simpleComputeTrans(G, rock);
    
    %% Fluid model
    
    % We set up a two-phase fluid. Viscosity, density is set for oil and water
    
    %gravity off
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

   s = linspace(0, 1, 1001).'; kr = fluid.relperm(s);
  f(2) = figure(2);
  plot(s, kr); legend('kr_1', 'kr_2'),title(['rel perm']), axis equal tight



    %Boundary conditions. We set boundary conditions, waterflooding

    pv = poreVolume(G, rock);
    %injRate = -sum(pv)/(500*day);
     injRate = -0.1*meter^3/day;
%     bc = fluxside([], G, 'ymin', -injRate, 'sat', [1, 0]);
%     bc = pside(bc, G, 'ymax', 0*barsa, 'sat', [0, 1]);
    
    W = addWell([],G, rock, 1:nx*ny:nx*ny*(nz), 'Type', 'bhp', ...
    'Val', 100*barsa, 'name', 'I', 'radius', .1, 'Comp_i', [1 0]);
   W = addWell(W,G, rock, nx*ny:nx*ny:G.cells.num, 'Type', 'bhp', ...
    'Val', 0, 'name', 'P', 'radius', .1, 'Comp_i', [0 1]);
%     
    % fluid = initSingleFluid('mu' , 1, ...
    %                         'rho', 1);
    % display(fluid);
%       plotGrid(G, 'FaceAlpha', 0, 'EdgeAlpha', .1); 
% figure
%       plotCellData(G,rock.perm)
%       hp=plotWell(G, W);
% %    
%     view(0,90)
%     plotGrid(G, 'FaceColor', 'none');

figure(1)

 plotWell(G, W(1),'color', 'r');
 plotWell(G, W(2), 'color', 'b');
            view(50,15) 
            colormap(jet) 
            axis equal tight 
   %return
%     
       %% Set pressure solver
        rSol = initState(G, [], 100*barsa, [0 1]);
        psolveropw

   
    %S  = computeMimeticIP(G, rock, 'Verbose', verbose,'InnerProduct','ip_tpf');
    fn = @(A, b) solver.solveLinearSystem(A, b);
    psolve = @(state) incompTPFA_g_o(state, G, hT, fluid, 'wells', W ,'MatrixOutput',true,'LinSolve', fn);
        %% Define transport solver
    
    switch tsolver
        case 1
            % Explicit tranport solver
            tsolve  = @(state, dT, fluid) explicitTransport(state, G, dT, rock, fluid,  'wells', W,  'verbose', false);
        case 2
            % Implicit transport solver: try with one time step
            tsolve  = @(state, dT, fluid) implicitTransport(state, G, dT, rock, fluid,  'wells', W, 'Verbose', false);
    end
    
    %% Initiate pressure solver


    heading = [' saturation'];
    %heading = [num2str(t),  ' days'];
    f(3) = figure(3);
    subplot('position',[0.001, 0, 0.5, 0.9]), cla
    
    plotCellData(G, rSol.s(:,1),'linestyle','none');
    caxis([0 1]), view(50,15), axis equal tight, title(['Water' heading]), colormap(jet),colorbar
    subplot('position',[0.501, 0, 0.5, 0.9]), cla
    plotCellData(G, rSol.s(:,2),'linestyle','none');
    caxis([0 1]), view(50,15), axis equal tight, title(['Oil' heading]), colormap(jet), colorbar
    
 %% Solve initial pressure in reservoir
    rSol = psolve(rSol);
         f(4) = figure(4);
    heading = [' pressure'];
    %heading = [num2str(t),  ' days'];
    
    subplot('position',[0.001, 0, 0.5, 0.9]), cla
    plotCellData(G, rSol.s(:,1),'linestyle','none');
    caxis([0 1]), view(50,15), axis equal tight off, title(['Water' heading]), colormap(jet),colorbar
    subplot('position',[0.501, 0, 0.5, 0.9]), cla
    plotCellData(G, rSol.s(:,2),'linestyle','none');
    caxis([0 1]), view(50,15), axis equal tight off, title(['Water' ' flux']), colormap(jet), colorbar
   %% Transport loop
   N      = fix(T/dTplot);
    pv     = poreVolume(G,rock);
    t      = 0;
    plotNo = 1;
    H1 = 'Water Sat - '; H2 = 'Oil Sat - '; H3 = 'W P - ';
    e = []; p_org = []; p_pc = [];
    figure;
    %dT=20;
    ts=0;
    podi=0;
    while t < T,
       
        ts=ts+1;
        % Increase time and continue if we do not want to plot saturations
%                   figure(1000)
%             pp=rSol.pressure;    
            if t==0
                clear p0
            p0=rSol.pressure;
            cmin=(4)*min(p0);
            cmax=0.7*max(p0);
            end
%         
%             plotGrid(G, 'FaceAlpha', 0, 'EdgeAlpha', .1);
%             plotCellData(G,pp)
%             axis equal tight off, colormap(jet), title('P0'),caxis([cmin cmax])
%             pause
%    
        
        %convertTo(t,day)
        % pause
        % TRANSPORT SOLVE
        [rSol,treport(ts)]   = tsolve(rSol, dT, fluid);
        %rSol_pc = tsolve(rSol_pc, dT, fluid_pc);
        
        % Check for inconsistent saturations
        %s = [rSol.s(:,1); rSol_pc.s(:,1)];
        s = [rSol.s(:,1)];
        assert(max(s) < 1+eps && min(s) > -eps);
        
        

        
        
        if  def==0
            [rSol,preport(ts)]    = psolve(rSol);
            
        else
            Zp(:,ts)=rSol.pressure;
            if ts <dv+1
                % Update solution of pressure equation. If backslash is used, there
                % is no report
                [rSol,preport(ts)]    = psolve(rSol);
                %rSol_pc = psolve(rSol_pc, fluid_pc);
                
            else

                podi=podi+1;
                Z=Zp(:,podi:podi+dv-1);
                
               % pause
                if pod==1
                    
                    [U,S]=defpodf_Dt(Z,dir2,dv,ts,dTplot/day());
                    Z=U(:,dpod);
                end
                nwcells = 0;
%                 for i = 1 : numel(W)
%                     nwcells = nwcells + numel(W(i).cells);
%                 end
                 for i=1:numel(W)
                            i
                            Z(nx*ny*nz+i,:)=0;
                  end
               lsolver = 4;
                 psolveropw
                fn = @(A, b) solver.solveLinearSystem(A, b);
                psolve = @(state) incompTPFA_g_o(state, G, hT, fluid,  'wells', W,'MatrixOutput',true,'LinSolve', fn);
                [rSol,preport(ts)]    = psolve(rSol);
                
            end
        end
       % condn(ts,1) = condest(rSol.A);
  
        t = t + dT;
 
        if ( t < plotNo*dTplot && t<T) continue, end
        % Plot saturation
       f(5) = figure(5);
        heading = [num2str(convertTo(t,day)),  ' days'];
        %heading = [num2str(t),  ' days'];
        r = 0.01;
        subplot('position',[(plotNo-1)/N+r, 0.50, 1/N-2*r, 0.4]), cla
        plotCellData(G, rSol.s(:,1),'linestyle','none');
        caxis([0 1]), view(50,15), axis equal tight off, title([H1 heading]), colormap(jet)
        subplot('position',[(plotNo-1)/N+r, 0.02, 1/N-2*r, 0.4]), cla
        plotCellData(G, rSol.s(:,2),'linestyle','none');
        caxis([0 1]), view(50,15), axis equal tight off, title([H2 heading]), colormap(jet)
        

        %heading = [num2str(t),  ' days'];
        r = 0.01;
        
        if t==dTplot;
          f(6) = figure(6);
        heading = [num2str(convertTo(dT,day)),  ' days'];
        subplot('position',[0.1, (plotNo)/(N+1)+r, 0.9, 1/(N+1)-r]), cla
        plotCellData(G, p0(1:n),'linestyle','none');
        view(50,15), %colorbar,
        axis equal tight off, colormap(jet), title([H3 heading],'Position',[(plotNo-1)/N+r+150 0.10 ]),caxis([cmin cmax])
        end
        %pause       
        f(6) = figure(6);
        heading = [num2str(convertTo(t,day)),  ' days'];
        p=rSol.pressure;
        subplot('position',[0.1, (plotNo)/(N+1)+r, 0.9, 1/(N+1)-r]), cla
        plotCellData(G, p,'linestyle','none');
        view(50,15), %colorbar,
        axis equal tight off, colormap(jet), title([H3 heading],'Position',[(plotNo-1)/N+r+150 0.10 ]),caxis([cmin cmax])
        plotNo = plotNo+1;
    end
    
   
%     
%     p=rSol.pressure;
% %     
%     A=rSol.A(1:G.cells.num,1:G.cells.num);
%     b=rSol.rhs(1:G.cells.num);
%     
%     % clf;
%     xb=A\b;
%     resb=p-xb;
    %[xd,flag,res,its]=DICCG_01_25_2(A,b,Z,tol,1000);
    %resd=xd-xb;
    
%     figure(per+1200)
%     title('Pressure')
%     [ht]=plotingsolution_bc(G,'(D)ICCG', p,1) ;
%     colorbar
%     [ht]=plotingsolution_bc(G,'backslash',xb,2);
%     colorbar
%     figure(per+1300)
%     plotGrid(G, 'FaceAlpha', 0, 'EdgeAlpha', .1);
%     plotCellData(G,xb-p)
%     view(0,90)
%     
%     axis equal tight; colormap(jet(128));
%     title('bs-DICGCG');
%     colorbar
    for j=1:ts
    its(j)=preport(j).iterations;
    end
    f(7) = figure(7);
    if def==1
        hn=plot(1:dv,its(1:dv),'r*');
        hold on
        hn=plot(dv+1:ts,its(dv+1:ts),'bp');
        hold on
        legend('ICCG','DICCG');
        axis square
    else
        hn=plot((dT:dT:T)/day,its(1:ts),'r*');
        legend('ICCG');
        axis square
    end
    ylabel('Number of iterations','FontSize',16)
    xlabel('Time (days)','FontSize',16)
    if def==0
    title(['Iterations ICCG'],'FontSize',16)
    else if (def==1)&&(pod==0)
    title(['Iterations DICCG, ' num2str(dv) ' deflation vectors'],'FontSize',16)
        else
    title(['Iterations DICCG, ' num2str(dv) ' snapshots, ' num2str(numel(dpod)) ' POD vectors '],'FontSize',16)    
        end
    end
%     f(8) = figure(8);
%             hn=plot((dT:dT:T)/day,condn(1:ts),'*');
%        % legend('ICCG');
%         axis square
%     ylabel(' \kappa_1(A)','FontSize',16)
%     xlabel('Time (days) ','FontSize',16)
%     if cn == 1
%      f(9) = figure(9);
%             hn=plot((dT:dT:T)/day,condn(1:ts),'*');
%        % legend('ICCG');
%         axis square
%     ylabel(' \kappa_1((P)M^{-1}A)','FontSize',16)
%     xlabel('Time (days) ','FontSize',16)
%     end
if cn == 1
defcondnumb
end
    %%If we want to save the files/graphs, it's neccesary to run savefilesf
    savefilesf

   
         
        end
    end
        
   end
end


