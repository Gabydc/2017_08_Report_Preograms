%%
% Two-phase Incompressible Flow Solver - Gabriela Diaz 17/05/17
% This program simulates two-phase flow through a porous media. The
% pressure equation is solved using DICCG Method and is compared with ICCG
close all
clear all
clc
% The capillary pressure can be taking into accoint (cp =1)
for cp = [0 ]
    
    % The contrast between permeability layers can be varied one layer has
    % permeability of 10*milli*darcy, the secon is varied as
    % 10*10^(-per)*milli*darcy
    for per=[0]
        
        % In this part the solver is chosen 0 means ICCG (no deflation) and
        % 1 is with deflation
        for def=[0 1 ]
            
            % If we want to use POD as deflation vectors pod=1
            for pod = [0]
                clearvars -except per def pod cp
                
                
                dpod=[];
                if (def == 0)  && (pod == 1)
                    pod=0;
                    continue
                end
                
                %We choose the number of POD deflation vectors to use. This number can be
                %varied in the vars subroutine
                for rr=1:2
                    if (pod == 0) && (rr==2)
                        continue
                    end
                    % Initialize the variables, size, fluid values, time
                    % step, ...
                    varswr
                    close all
                    dpod=podv{rr};
                    if pod== 0
                        dpod=[];
                    end
                    
                    
                    
                    %Create the directory
                    dir='/mnt/sda2/cortes/Results/17_09/training/';
                    
                    folder=[ '10-' num2str(k) '_' num2str(sz) 'nz' num2str(nz) 'perm_' num2str(per) 'cp' num2str(cp)];
                    mkdir([dir], folder)
                    dir1 = [dir folder '/'];
                    
                    folder=[  'def_' num2str(def) '_pod_' num2str(numel(dpod)) ];
                    mkdir([dir1], folder)
                    dir2 = [dir1 folder '/'];
                    
                    %Create the grid
                    G = cartGrid([sz, sz, nz], [Lx, Ly, Lz]);
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
                 %     well(1:4)=4;
                   well(5)=795;  
                    well(5) = 0.3;
wtype    = {'bhp', 'bhp', 'bhp', 'bhp', 'bhp'};
wtarget  = [well(1),   well(2),   well(3),   well(4), well(5)/(barsa()*day())] .* barsa();
wrad     = [0.125, 0.125, 0.125, 0.125, 0.125] .* meter;
wloc     = [  nxi,   nxi,     nx,   nx, nx/2;
              nyi,   ny,     nyi,   ny, ny/2];
%wname    = {'I1', 'I2', 'I3', 'I4', 'P'};
wname    = {'P1', 'P2', 'P3', 'P4', 'I'};
sgn      = [ -1 ,  -1 ,  -1 ,  -1 , 1 ];
Comp_iI= [1 0];
Comp_iP= [0 1];

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
                    
                    % Create a well structure to support multiple report steps
                    [newW{1:steps}] = deal(W);
                    BHP = zeros(numel(W),steps);
                    Valr = zeros(numel(W),steps);
                    %Then we'll assign the time-dependent controls
                    BHP(5,:) = linspace(20, 80, steps)*barsa;
%                     Valr(1,1:steps/2) = linspace(0.5, 0.8, steps/2)/day;
%                     Valr(2,steps/2+1:steps) = linspace(0.5, 0.8, steps/2)/day;
%                     Valr(4,1:steps/2) = linspace(0.5, 0.8, steps/2)/day;
%                     Valr(3,steps/2+1:steps) = linspace(0.5, 0.8, steps/2)/day;
                    Valr(5,1:steps) = linspace(0.5, 0.8, steps)/day;
                    %Valr =  linspace(0.5, 1, steps)/day();
                    for w = 1:numel(newW),
                        for i = 5
                          %  newW{w}(i).type = 'bhp';
                          %  newW{w}(i).val = BHP(i,w);
                            newW{w}(i).type = 'rate';
                           newW{w}(i).val = Valr(i,w);
                        end
                    end
                    
                    
                    f(1)=figure(1);
                    plotCellData(G, rock.perm/(milli*darcy())); colorbar
                    colormap(jet), axis equal tight off
for i = 1 :numel(W)-1
                    plotWell(G, W(i),'color', 'r');
end
                    plotWell(G, W(5), 'color', 'b');
                    title('Permeability field [mD]')
                    if nz > 1
                        view(10,20)
                    else
                        view(0,90)
                    end
                    
                   %return
                    %% Set pressure solver
                    % Initiate the pressure of the reservoir
                    rSol = initState(G, [], 0*barsa, [0 1]);
                    IWS =  rSol.s(:,1);
                    IWP = rSol.pressure;
                    % Set linear solver for incompTPFA_g_o ( The difference with incompTPFA is
                    % that we can get the report
                    psolveropw
                    fn = @(A, b) solver.solveLinearSystem(A, b);
                    psolve = @(state) incompTPFA_g_o(state, G, hT, fluid, 'wells', W,'MatrixOutput',true,'LinSolve', fn);
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
                    % Solve initial pressure in reservoir
                    rSol = psolve(rSol);
                    
                    %% Transport loop
                    N      = fix(T/dTplot);
                    pv     = poreVolume(G,rock);
                    t      = 0;
                    plotNo = 1;
                    H1 = 'W Sat - '; H2 = 'Oil Sat - '; H3 = 'Pressure ';
                    e = []; p_org = []; p_pc = [];
                    ts=0;
                    podi=0;
                    while t < T,
                        ts=ts+1;
                        if t==0
                            p=rSol.pressure;
                            cmin=(4)*min(p);
                            cmax=0.7*max(p);
                        end
                                                       % Pressure Solver
                                
                                W = newW{ts};
                                fn = @(A, b) solver.solveLinearSystem(A, b);
                                psolve = @(state) incompTPFA_g_o(state, G, hT, fluid,  'wells', W,'MatrixOutput',true,'LinSolve', fn);
                               
                        % TRANSPORT SOLVE
                        [rSol,treport(ts)]   = tsolve(rSol, dT, fluid);
                        % Check for inconsistent saturations
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
                            else
                                podi=podi+1;
                                Z=Zp(:,podi:podi+dv-1);
                                if pod==1
                                    
                                    [U,S]=defpodf_Dt(Z,dir2,dv,t/day(),dTplot/day());
                                    Z=U(:,dpod);
                                end
                                nwcells = 0;
                                for i = 1 : numel(W)
                                    nwcells = nwcells + numel(W(i).cells);
                                end
                                for i=1:numel(W)
                                    
                                    Z(nx*ny*nz+i,:)=0;
                                end
                                lsolver = 4;
                                psolveropw
                                W = newW{ts};
                                fn = @(A, b) solver.solveLinearSystem(A, b);
                                psolve = @(state) incompTPFA_g_o(state, G, hT, fluid,  'wells', W,'MatrixOutput',true,'LinSolve', fn);
                                [rSol,preport(ts)]    = psolve(rSol);
                                
                            end
                        end
                        if ceigs == 1
                            [Vn,Dn] = eigs(rSol.A,length(rSol.A));
                            Dn = diag(Dn);
                            Dn = abs(Dn);
                            digits(4)
                            condn(ts,1)=max(Dn)/min(Dn);
                        end
                        wellSols{ts,1}=getWellSol(W, rSol, fluid);
                        t = t + dT;
                        
                        if ( t < plotNo*dTplot && t<T) continue, end
                        
                        % Plot saturation
                        f(3) = figure(3);
                        heading = [num2str(convertTo(t,day)),  ' days'];
                        r = 0.01;
                        % Plot intitial conditions
                        if plotNo == 1;
                            
                            subplot('position',[(plotNo-1)/(N+2)+r, 0.1, 1/(N+2)-r, 0.4]), cla
                            plotCellData(G, IWP,'LineStyle',':');
                            %caxis([cmin cmax]),
                             axis equal tight off, title([H3 ]), colormap(jet)
                            if nz > 1
                                view(-25,20)
                            else
                                view(0,90)
                            end
                            subplot('position',[(plotNo-1)/(N+2)+r, 0.5, 1/(N+2)-r, 0.4]), cla
                            plotCellData(G, IWS,'LineStyle',':');
                            caxis([0 1]), axis equal tight off, title([H1 '0 days']), colormap(jet)
                            if nz > 1
                                view(-25,20)
                            else
                                view(0,90)
                            end
                        end
                        % Plot rest of the time steps
                        subplot('position',[(plotNo)/(N+2)+r, 0.5, 1/(N+2)-r, 0.4]), cla
                        plotCellData(G, rSol.s(:,1),'LineStyle',':');
                        caxis([0 1]), axis equal tight off, title([ heading]), colormap(jet)
                        if nz > 1
                            view(-25,20)
                        else
                            view(0,90)
                        end
                        if t==T;
                            caxis([0 1]), axis off
                            colorbar('Position', [(plotNo+1)/(N+2)+0.05, 0.6, 4*r, 0.2])
                            if nz > 1
                                view(-25,20)
                            else
                                view(0,90)
                            end
                        end
                        %   (plotNo-1)/N+r, 0.02, 1/N-2*r, 0.4]
                        subplot('position',[(plotNo)/(N+2)+r, 0.1, 1/(N+2)-r, 0.4]), cla
                        plotCellData(G, rSol.pressure(:),'LineStyle',':');
                       % caxis([cmin cmax]),
                          axis equal tight  off,  colormap(jet)
                        if nz > 1
                            view(-25,20)
                        else
                            view(0,90)
                        end
                        
                        
                        if t==T;
                            %caxis([cmin cmax]), 
axis off
                            colorbar('Position', [(plotNo+1)/(N+2)+0.05, 0.2, 4*r, 0.2])
                            if nz > 1
                                view(-25,20)
                            else
                                view(0,90)
                            end
                        end
                        plotNo = plotNo+1;
                    end
                    
                    for j=1:ts
                        its(j)=preport(j).iterations;
                    end
                    f(4) = figure(4);
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
                    
                    if cn == 1
                        defcondnumb
                    end
                    %%If we want to save the files/graphs, it's neccesary to run savefilesf
                    % savefilesf
                    
                    
                    
                end
            end
        end
    end
end

                nf = nf + 1;
                figure(nf);
                file{nf} = ['Production_curves'];
                plot(t,cellfun(@(x) x(2).Sw, wellSols), ...
                    t,cellfun(@(x) x(2).wcut, wellSols));
                legend('Sw in completion','Water cut','Location','NorthWest');
                axis([0 max(t) -.05 1.0]);
                xlabel('Time (PVI)','FontSize',14)
                ylabel('Saturation','FontSize',14)
                
                %%
                % Second, we plot the oil rate used in our simulation in units m^3/day
                nf = nf + 1;
                figure(nf);
                file{nf} = ['Oil_rate'];
                
                for i = 1 : numel(W)
                    qOs = cellfun(@(x) abs(x(i).qOs), wellSols);
                    plot(t, qOs([2:end end])*day,'LineStyle','--', 'Color', pcol{i},'Marker',pmark(i),'MarkerSize',4); %axis([0 1.2 ]);
                    hold on
                end
                legend(wname,'Location','NorthEast');
                axis tight
                xlabel('Time (PVI)','FontSize',14)
                ylabel('Oil rate [m^3/day]','FontSize',14)
                %%
                % Second, we plot the water rate used in our simulation in units m^3/day
                nf = nf + 1;
                figure(nf);
                file{nf} = ['Water_rate'];
                for i = 1 : numel(W)
                    qWs = cellfun(@(x) abs(x(i).qWs), wellSols);
                    plot(t, qWs([2:end end])*day,'LineStyle','--', 'Color', pcol{i},'Marker',pmark(i),'MarkerSize',4); %axis([0 1.2 ]);
                    hold on
                end
                legend(wname,'Location','NorthEast'); axis tight
                xlabel('Time (PVI)','FontSize',14)
                ylabel('Water rate [m^3/day]','FontSize',14)
                %
                % Second, we plot the water rate used in our simulation in units m^3/day
                nf = nf + 1;
                figure(nf);
                file{nf} = ['bhp'];
                for i = 1 : numel(W)
                    bhp = cellfun(@(x) abs(x(i).bhp), wellSols);
                    bhpv = bhp([2:end end]);
                    plot(t, bhpv./barsa,'LineStyle','--', 'Color', pcol{i},'Marker',pmark(i),'MarkerSize',4); %axis([0 1.2 ]);
                    hold on
                end
                legend(wname,'Location','NorthEast'); axis tight
                xlabel('Time (PVI)','FontSize',14)
                ylabel('Bhp [bars]','FontSize',14)
                
                
                %
                % Last, we plot the cumulative oil production computed from the well
                % solution and compare this with the amount of extracted oil derived from a
                % mass-balance computation (initial oil in place minus current oil in
                % place). We also include a horizontal line indicating the initial oil in
                % place, and a straight line showing the amount of oil we would have
                % extracted if oil was produced at the constant initial rate
%                 nf = nf + 1;
%                 figure(nf);
%                 file{nf} = ['Cum_Oil_P'];
%                 plot(t,cumsum(bsxfun(@times, abs(cellfun(@(x) x(2).qOs, wellSols)), dt)));
%                 hold on;
%                 plot(t,oip(1)-oip,'-.','LineWidth',3);
%                 plot([0 1.2],oip([1 1]),'-k',t,min(t*oip(1),oip(1)),'-k', ...
%                     t(W(2).bt+[0 0]),[0 oip(1)],'--k');
%                 hold off; axis tight; axis([0 max(t) 0 1.05*oip(1)]);
%                 xlabel('Time (PVI)','FontSize',14)
%                 ylabel('Cumulative Oil Production [m^3/day]','FontSize',14)
                
                
                % %% Plotting with GUI from ad-core
                % mrstModule add ad-core
                % plotWellSols(wellSols,cumsum(dt))
                %%
                if sf == true
                    for i = 1 : nf
                        f(i) = figure(i);
                        savefigures(f(i), file{i}, dir2)
                    end
                end
