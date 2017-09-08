close all
clear all
clc





% for i=1:7
%     f(i) = figure('visible','off');
% end
for cp = [0 ]
    cp
    %training = 0, training phase
    %training = 1, test cases
    training =[1 ]
    for def=[0  1]
        def

        % If we want to use POD pod==1
        
      %  for pod = [1]
            pod = 1;
            dpod=[];
            
            
%             if (def == 0)  && (pod == 1)
%                 dpod = [];
%                 pod = 0;
%                 continue
%             end
%             
            
            
            for rr=1:2
                rr
                
                if (pod == 0) && (rr == 2)
                    continue
                end
                if (pod == 0) && (def == 1)
                    continue
                end
                close all
                clearvars -except per def pod rr dpod cp training
                variablesSPE10_6_t
                
                
                
                
                
                
              
                
                %% Compute half transmissibilities
                % All we need to know to develop the spatial discretization is the reservoir
                % geometry and the petrophysical properties. This means that we can compute
                % the half transmissibilities without knowing any details about the fluid
                % properties and the boundary conditions and/or sources/sinks that will
                % drive the global flow:
                hT = computeTrans(G, rock);
                
                
                
                
                %% Initiate pressure solver
                rSol = initState(G, W, 500*barsa, [0 1]);
                
                %% Set pressure solver
                psolverop2w
                
                %% Solve initial pressure in reservoir
                %rSol = psolve(rSol);
                rSol = incompTPFA(rSol, G, hT, fluid, 'wells', W);
                IWS =  rSol.s(:,1);
                IWP = rSol.pressure;
                
                %% Transport loop
                %  T      = 100*day();
                T  = 1.2*sum(pv)/rSol.wellSol(1).flux;
                % T = max(T);
                % T = 10000*day;
                steps=480;
                %steps =48;
                Stot = steps;
                dT     = ceil(T/steps);
                dTplot = T/4;  % plot only every 100th day
                N      = fix(T/dTplot);
                per =0;
                t      = 0;
                plotNo = 1;
                H1 = 'Water Sat - '; H2 = 'Oil Sat - '; H3 = 'W P - ';
                e = []; p_org = []; p_pc = [];
                %figure;
                %dT=20;
                ts=0;
                podi=0;
                
                
                % Compute time-of-flight for the single-phase flow field and record the
                % corresponding breakthrough time in the producer.
                tau = computeTimeOfFlight(rSol, G, rock,  'wells', W);
                tbf = tau(W(2).cells,1);
                
                % Initialize number of time intervals, cell array to hold well solutions,
                % and array to hold the oil in place
                wellSols = cell(N+1,1);  wellSols{1} = getWellSol(W, rSol, fluid);
                oip      = zeros(N,1); oip(1) = sum(rSol.s(:,2).*pv);
                while t < T,
                    
                    ts=ts+1
                    if training == 0
                    W =W0{ts,1};
                    else
                    W =W1{ts,1};
                    end
                    
                    % Increase time and continue if we do not want to plot saturations
                    
                    if t==0
                        p0=rSol.pressure;
                        
                    end
                    p0 = rSol.pressure;
                    n=nx*ny*nz;
                    for i=1:numel(W)
                        p0(n+i) = 0;
                    end
                    if training == 1
                        load([dir1 'Z.mat'],'Z')  
                        
                    end
                    if  def==0
                        solver = PCG_ICSolverAD_cn('tolerance', tol,'maxIterations', maxIterations,'cn',cn);
                        fn = @(A, b) solver.solveLinearSystem(A, b);
                        psolve = @(state) incompTPFA_g_o(state, G, hT, fluid, 'wells', W, 'MatrixOutput',true,'LinSolve', fn);
                        [rSol,preport(ts)]    = psolve(rSol);
                        
                    else                                    
                            solver = DPCG_ICSolverAD('tolerance', tol,'maxIterations', maxIterations, 'Z',Z,'x0',p0);
                            fn = @(A, b) solver.solveLinearSystem(A, b);
                            psolve = @(state) incompTPFA_g_o(state, G, hT, fluid, 'wells', W, 'MatrixOutput',true,'LinSolve', fn);
                            [rSol,preport(ts)]    = psolve(rSol);
                        end
                        if ceigs == 1
                            [Vn,Dn] = eigs(rSol.A,length(rSol.A));
                            Dn = diag(Dn);
                            Dn = abs(Dn);
                            digits(4)
                            condn(ts,1)=max(Dn)/min(Dn);
                        end
                        
                        
                        % TRANSPORT SOLVE
                        [rSol]   = tsolve(rSol, dT, fluid);
                        %rSol_pc = tsolve(rSol_pc, dT, fluid_pc);
                        wellSols{ts+1} = getWellSol(W, rSol, fluid);
                        oip(ts+1) = sum(rSol.s(:,2).*pv);
                        
                        if rSol.s(W(2).cells,1)<eps
                            W(2).bt = ts;
                        end
                        Sat_w(:,ts) =  rSol.s(:,1);
                        Pressure(:,ts) = rSol.pressure(:)/barsa;
                        Pressure1(:,ts) = rSol.pressure(:)/barsa;
                        % Check for inconsistent saturations
                        %s = [rSol.s(:,1); rSol_pc.s(:,1)];
                        s = [rSol.s(:,1)];
                        assert(max(s) < 1+eps && min(s) > -eps);
                        
                        
                        t = t + dT;
                        
                        
                        
                        
                   
                    
                end
                
                if training == 0
                    
                    Z =  Pressure1;
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
                end
                
              
                
                
               


                pmin=min(Pressure1(:,1));
                pmax=max(Pressure1(:,ts));
                
                plotNo = 0;
                Np = 3;
                time = (0:dT:T)/day;
                px = [0.01 0.26 0.51 0.76];
                
                nf = nf + 1;
                figure(nf);
                %title('Water Saturation');
                file{nf} = ['Water_saturation'];
                clim = [0 1];
                subplotcbspe(nf,clim,ts,Np,G,nz,time,Sat_w)
                % subplotcbwbt(nf,clim,wbt,Np,G,nz,time,Sat_w)
                nf1 = nf + 1;
                figure(nf1);
                %title('Pressure [bars]');
                file{nf1} = ['Pressure'];
                clim = [pmin pmax];
                subplotcbspe(nf1,clim,ts,Np,G,nz,time,Pressure1)
                % subplotcbwbt(nf1,clim,wbt,Np,G,nz,time,Pressure1)
                nf=nf1;
                
                %
                for j=1:ts
                    its(j)=preport(j).iterations;
                end
                nf = nf + 1;
                f(nf) = figure(nf);
                file{nf} = ['Iterations'];
                if def==0
                    plot(time(1:ts),its(1:ts),'r*');
                    %title(['Iterations ICCG'],'FontSize',16)
                    legend('ICCG');
                    axis([ dT/day T/day 0 max(its)]);
                    ylabel('Number of iterations','FontSize',16)
                    xlabel('Time (days)','FontSize',16)
                    axis square
                else if (def==1)&&(pod==0)
   
                        plot(time(1:ts),its(1:ts),'bp');
                        legend('DICCG');
                        axis([ dT/day T/day 0 max(its)]);
                        ylabel('Number of iterations','FontSize',16)
                        xlabel('Time (days)','FontSize',16)
                        axis square
                        % title(['Iterations DICCG, ' num2str(dv) ' deflation vectors'],'FontSize',16)
                    else
 
                        plot(time(1:ts),its(1:ts),'bp');
                        legend(['DICCG_{POD_{' num2str(numel(dpod)) '}}' ]);
                        axis([ dT/day T/day 0 max(its)]);
                        ylabel('Number of iterations','FontSize',16)
                        xlabel('Time [days]','FontSize',16)
                        axis square
                        % title(['Iterations DICCG, ' num2str(dv) ' snapshots, ' num2str(numel(dpod)) ' POD vectors '],'FontSize',16)
                    end
                end
                
                
                
                
                
if training == [1 ]
filetx = ['results.txt'];
                saveres(dir1,filetx,def,pod,dpod,per,ts,dv,preport)
end

                
                savefilesfun1
                %     %% Plot production curves
                % % First we plot the saturation in the perforated cell and the corresponding
                % % water cut, i.e., the fractional flow evaluated in the completion
                dt = [0; ones(ts,1)*T/(ts)]; t = 1.2*cumsum(dt)/T;
                
                
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
                
                
                                nf = nf + 1;
                figure(nf);
                file{nf} = ['bhp_p'];
                for i = 2 : numel(W)-1
                    bhp = cellfun(@(x) abs(x(i).bhp), wellSols);
                    bhpv = bhp([2:end end]);
                    plot(time(1:ts), bhpv(1:ts)./barsa,'LineStyle','--', 'Color', pcol{i},'Marker',pmark(i),'MarkerSize',4); %axis([0 1.2 ]);
                    hold on
                end
                legend(wname(2:5),'Location','NorthEast'); axis tight
                xlabel('Time [days]','FontSize',14)
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
                %
                
            end
        %end
        
    end
    
    end









