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
        for def=[0 ]
            
            % If we want to use POD as deflation vectors pod=1
            for pod = [0 1]
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
                    vars
                    close all
                    dpod=podv{rr};
                    if pod== 0
                        dpod=[];
                    end
                    
                    
                    
                    %Create the directory
                    dir='/mnt/sda2/cortes/Results/2017/Report/bc/';
                    
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
                    f(1)=figure(1);
                    plotCellData(G, rock.perm/(milli*darcy())); colorbar
                    colormap(jet), axis equal tight off
                    title('Permeability field [mD]')
                    if nz > 1
                        view(10,20)
                    else
                        view(0,90)
                    end
                    %                     break
                    
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
                    injRate = -0.4*meter^3/day;
                    bc = fluxside([], G, 'xmin', -injRate, 'sat', [1, 0]);
                    bc = pside(bc, G, 'xmax', 0*barsa, 'sat', [0, 1]);
                    
                    %% Set pressure solver
                    % Initiate the pressure of the reservoir
                    rSol = initState(G, [], 100*barsa, [0 1]);
                    IWS =  rSol.s(:,1);
                    IWP = rSol.pressure;
                    % Set linear solver for incompTPFA_g_o ( The difference with incompTPFA is
                    % that we can get the report
                    psolverop
                    fn = @(A, b) solver.solveLinearSystem(A, b);
                    psolve = @(state) incompTPFA_g_o(state, G, hT, fluid, 'bc', bc ,'MatrixOutput',true,'LinSolve', fn);
                    %% Define transport solver
                    switch tsolver
                        case 1
                            % Explicit tranport solver
                            tsolve  = @(state, dT, fluid) explicitTransport(state, G, dT, rock, fluid,  'bc', bc,  'verbose', false);
                        case 2
                            % Implicit transport solver: try with one time step
                            tsolve  = @(state, dT, fluid) implicitTransport(state, G, dT, rock, fluid,  'bc', bc, 'Verbose', false);
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
                                lsolver = 4;
                                psolverop
                                fn = @(A, b) solver.solveLinearSystem(A, b);
                                psolve = @(state) incompTPFA_g_o(state, G, hT, fluid,  'bc', bc,'MatrixOutput',true,'LinSolve', fn);
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
                        t = t + dT;
                        
                        if ( t < plotNo*dTplot && t<T) continue, end
                        % Plot saturation
                        f(3) = figure(3);
                        heading = [num2str(convertTo(t,day)),  ' days'];
                        r = 0.01;
                        % Plot intitial conditions
                        if t==dTplot;
                            subplot('position',[(plotNo-1)/(N+2)+r, 0.1, 1/(N+2)-r, 0.4]), cla
                            plotCellData(G, IWP,'LineStyle','none');
                            caxis([cmin cmax]),axis equal tight off, title([H3 ]), colormap(jet)
                            if nz > 1
                                view(-25,20)
                            else
                                view(0,90)
                            end
                            subplot('position',[(plotNo-1)/(N+2)+r, 0.5, 1/(N+2)-r, 0.4]), cla
                            plotCellData(G, IWS,'LineStyle','none');
                            caxis([0 1]), axis equal tight off, title([H1 '0 days']), colormap(jet)
                            if nz > 1
                                view(-25,20)
                            else
                                view(0,90)
                            end
                        end
                        % Plot rest of the time steps
                        subplot('position',[(plotNo)/(N+2)+r, 0.5, 1/(N+2)-r, 0.4]), cla
                        plotCellData(G, rSol.s(:,1),'LineStyle','none');
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
                        plotCellData(G, rSol.pressure(:),'LineStyle','none');
                        caxis([cmin cmax]),axis equal tight  off,  colormap(jet)
                        if nz > 1
                            view(-25,20)
                        else
                            view(0,90)
                        end
                        
                        
                        if t==T;
                            caxis([cmin cmax]), axis off
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
                    savefilesf
                    
                    
                    
                end
            end
        end
    end
end


