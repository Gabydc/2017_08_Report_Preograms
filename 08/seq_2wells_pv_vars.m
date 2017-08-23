close all
clear all
clc





% for i=1:7
%     f(i) = figure('visible','off');
% end
for cp = [0 1]
    cp
for per=[0 1 2 3 4 5 6 7 ]
    per
for def=[0 1]
    def
    
    % If we want to use POD pod==1
    
    for pod = [0 1]

dpod=[];

          
        if (def == 0)  && (pod == 1)
            dpod = [];
            pod = 0;
            continue
        end

        
        
        for rr=1:2
            rr
            
            if (pod == 0) && (rr == 2)
                continue
            end 
            close all
        clearvars -except per def pod rr dpod cp
        variables2w
            

         
            
              
                
                clear Z
                
                %% Compute half transmissibilities
                % All we need to know to develop the spatial discretization is the reservoir
                % geometry and the petrophysical properties. This means that we can compute
                % the half transmissibilities without knowing any details about the fluid
                % properties and the boundary conditions and/or sources/sinks that will
                % drive the global flow:
                hT = computeTrans(G, rock);
                
                
                
                
                %% Initiate pressure solver
                rSol = initState(G, W, 100*barsa, [0 1]);
                
                %% Set pressure solver
                psolverop2w
                
                %% Solve initial pressure in reservoir
                %rSol = psolve(rSol);
                rSol = incompTPFA(rSol, G, hT, fluid, 'wells', W);
                IWS =  rSol.s(:,1);
                IWP = rSol.pressure;
                
                %% Transport loop
                %T      = 900*day();
                if per == 0
                T  = 1.2*sum(pv)/rSol.wellSol(1).flux;
                save T T
                else
                    load T
                end
                dT     = ceil(T/240);
                dTplot = T/4;  % plot only every 100th day
                N      = fix(T/dTplot);
                
                t      = 0;
                plotNo = 1;
                H1 = 'Water Sat - '; H2 = 'Oil Sat - '; H3 = 'W P - ';
                e = []; p_org = []; p_pc = [];
                %figure;
                %dT=20;
                ts=0;
                podi=0;
                nf = nf + 1;
                 figure(nf);
                title('Water Saturation');
                file{nf} = ['Water_saturation'];
                nf1 = nf + 1;
                 figure(nf1);
                title('Pressure [bars]');
                file{nf1} = ['Pressure'];
                
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
                    % Increase time and continue if we do not want to plot saturations
                    
                    if t==0
                        p0=rSol.pressure;

                    end
                    
                    
                    
                    
                  % rSol= incompTPFA(rSol, G, hT, fluid, 'wells', W);
                    
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
                            if (def == 1) && (pod == 0)
                                Z= eye(size(Z));
                            end
                            
                            if pod==1
                                
                                [U,S]=defpodf_Dt(Z,dir2,dv,ts,dTplot/day());
                                Z=U(:,dpod);
                            end
                            
                            for i=1:numel(W)
                                
                                Z(nx*ny+i,:)=0;
                            end
                            
                            solver = DPCG_ICSolverAD('tolerance', tol,'maxIterations', maxIterations, 'Z',Z);
                            fn = @(A, b) solver.solveLinearSystem(A, b);
                            psolve = @(state) incompTPFA_g_o(state, G, hT, fluid, 'wells', W, 'MatrixOutput',true,'LinSolve', fn);
                            [rSol,preport(ts)]    = psolve(rSol);
                        end
                    end
%                     

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
                
plotNo = 0;
Np = 4;
time = (dT:dT:T)/day;
px = [0.01 0.26 0.51 0.76];

for i = 1 : ceil(ts/Np) : ts
    i
                       % ptpv=t*1.2/T;
                       ptpv = 1.2*((plotNo)/(Np-1));
                       if per ==0
                     heading = [sprintf('t=%.2f PVI',ptpv)];
                       else
                           
                           heading = [sprintf('t=%d days',ceil(time(i)))];
                       end
                     figure(nf);
                    file{nf} = ['Saturation'];
                    
                    h=subplot(1,4,plotNo+1,'Position',[px(plotNo+1),0,0.2,0.8]);
                    axis equal tight off,  title([ heading])
                    plotCellData(G, Sat_w(:,i),'LineStyle','none');  
                    
                    
              
                    if nz > 1
                        view(-25,20)
                    else
                        view(0,90)
                    end
                    figure(nf1);
                    file{nf1} = ['Pressure'];
                    h1=subplot(1,4,plotNo+1,'Position',[px(plotNo+1),0,0.2,0.8]);
                    plotCellData(G, Pressure1(:,i),'LineStyle','none');
                     axis equal tight off, %caxis([cmin cmax])
                     title([ heading])
                     
                    if nz > 1
                        view(-25,20)
                    else
                        view(0,90)
                    end
                    plotNo = plotNo + 1;
end
nf=nf1;
%                    
                for j=1:ts
                    its(j)=preport(j).iterations;
                end
                nf = nf + 1;
                figure(nf);
                file{nf} = ['Iterations'];
                if def==1
                    
                    h=plot(1:dv,its(1:dv),'r*');
                    hold on
                    h=plot(dv+1:ts,its(dv+1:ts),'bp');
                    hold on
                    legend('ICCG','DICCG');
                    axis square
                else
                    h=plot(1:ts,its(1:ts),'r*');
                    legend('ICCG');
                    axis square
                end
                ylabel('Number of iterations','FontSize',16)
                xlabel('Time step ','FontSize',16)
                if def==0
                    title(['Iterations'],'FontSize',16)
                else if (def==1)&&(pod==0)
                        title(['Iterations, ' num2str(dv) ' deflation vectors'],'FontSize',16)
                    else
                        title(['Iterations, ' num2str(dv) ' snapshots, ' num2str(numel(dpod)) ' POD vectors '],'FontSize',16)
                    end
                end
                
                savefilesfun 
    %% Plot production curves
% First we plot the saturation in the perforated cell and the corresponding
% water cut, i.e., the fractional flow evaluated in the completion
dt = [0; ones(ts,1)*T/(ts)]; t = 1.2*cumsum(dt)/T; 
time = (0:dT:T+dT)/day;

nf = nf + 1;
figure(nf);
file{nf} = ['Production_curves'];
if per == 0
    plot(t,cellfun(@(x) x(2).Sw, wellSols),'--', ...
    t,cellfun(@(x) x(2).wcut, wellSols),'-','LineWidth',1);
legend('Sw in completion','Water cut','Location','NorthWest');
axis([0 max(t) -.05 1.0]);
 xlabel('Time [PVI]','FontSize',14)
else
      plot(time,cellfun(@(x) x(2).Sw, wellSols),'--', ...
    time,cellfun(@(x) x(2).wcut, wellSols),'-','LineWidth',1);
legend('Sw in completion','Water cut','Location','NorthWest');
axis([0 max(time) -.05 1.0]);
 xlabel('Time [days]','FontSize',14)  
end

ylabel('Saturation','FontSize',14)
  
%%
% Second, we plot the oil rate used in our simulation in units m^3/day
nf = nf + 1;
figure(nf);
file{nf} = ['Oil_rate'];
if per == 0
qOs = cellfun(@(x) abs(x(2).qOs), wellSols);
qOs1 = cellfun(@(x) abs(x(1).qOs), wellSols);
plot(t, qOs([2:end end])*day,'LineWidth',1); %axis([0 1.2 ]); 
hold on
plot(t, qOs1([2:end end])*day,'LineWidth',1);% axis([0 1.2 ]);
legend('P','I','Location','NorthEast');
axis tight
 xlabel('Time [PVI]','FontSize',14)
else
 qOs = cellfun(@(x) abs(x(2).qOs), wellSols);
qOs1 = cellfun(@(x) abs(x(1).qOs), wellSols);
plot(time, qOs([2:end end])*day,'LineWidth',1); %axis([0 1.2 ]); 
hold on
plot(time, qOs1([2:end end])*day,'LineWidth',1);% axis([0 1.2 ]);
legend('P','I','Location','NorthEast');
axis tight
 xlabel('Time [days]','FontSize',14)
end
ylabel('Oil rate [m^3/day]','FontSize',14)  
%%
% Second, we plot the water rate used in our simulation in units m^3/day
nf = nf + 1;
figure(nf);
file{nf} = ['Water_rate'];
if per == 0
qWs = cellfun(@(x) abs(x(2).qWs), wellSols);
qWs1 = cellfun(@(x) abs(x(1).qWs), wellSols);
plot(t, qWs([2:end end])*day,'LineWidth',1); %axis([0 1.2 ]);
hold on
plot(t, qWs1([2:end end])*day,'LineWidth',1);% axis([0 1.2 ]);
hold on
legend('P','I','Location','NorthEast'); axis tight
 xlabel('Time [PVI]','FontSize',14)
else
 qWs = cellfun(@(x) abs(x(2).qWs), wellSols);
qWs1 = cellfun(@(x) abs(x(1).qWs), wellSols);
plot(time, qWs([2:end end])*day,'LineWidth',1); %axis([0 1.2 ]);
hold on
plot(time, qWs1([2:end end])*day,'LineWidth',1);% axis([0 1.2 ]);
hold on
legend('P','I','Location','NorthEast'); axis tight
 xlabel('Time ([days]','FontSize',14)
end
ylabel('Water rate [m^3/day]','FontSize',14)
%%
% Second, we plot the water rate used in our simulation in units m^3/day
nf = nf + 1;
figure(nf);
file{nf} = ['bhp'];
if per == 0
bhp = cellfun(@(x) abs(x(2).bhp), wellSols);
bhp1 = cellfun(@(x) abs(x(1).bhp), wellSols);
bhpv = bhp([2:end end]);
bhpv1 = bhp1([2:end end]);
plot(t, bhpv./barsa,'LineWidth',1); %axis([0 1.2 ]);
hold on
plot(t, bhpv1./barsa,'LineWidth',1);% axis([0 1.2 ]);
legend('P','I','Location','NorthEast'); axis tight
 xlabel('Time [PVI]','FontSize',14)
else
  bhp = cellfun(@(x) abs(x(2).bhp), wellSols);
bhp1 = cellfun(@(x) abs(x(1).bhp), wellSols);
bhpv = bhp([2:end end]);
bhpv1 = bhp1([2:end end]);
plot(time, bhpv./barsa,'LineWidth',1); %axis([0 1.2 ]);
hold on
plot(time, bhpv1./barsa,'LineWidth',1);% axis([0 1.2 ]);
legend('P','I','Location','NorthEast'); axis tight
 xlabel('Time [days]','FontSize',14) 
end
ylabel('Bhp [bars]','FontSize',14)  


%%
% Last, we plot the cumulative oil production computed from the well
% solution and compare this with the amount of extracted oil derived from a
% mass-balance computation (initial oil in place minus current oil in
% place). We also include a horizontal line indicating the initial oil in
% place, and a straight line showing the amount of oil we would have
% extracted if oil was produced at the constant initial rate
nf = nf + 1;
figure(nf);
file{nf} = ['Cum_Oil_P'];
if per ==0
plot(t,cumsum(bsxfun(@times, abs(cellfun(@(x) x(2).qOs, wellSols)), dt)));
hold on;
plot(t,oip(1)-oip,'-.','LineWidth',3);
plot([0 1.2],oip([1 1]),'-k',t,min(t*oip(1),oip(1)),'-k', ...
    t(W(2).bt+[0 0]),[0 oip(1)],'--k');
hold off; axis tight; axis([0 max(t) 0 1.05*oip(1)]);
 xlabel('Time [PVI]','FontSize',14)
else
plot(time,cumsum(bsxfun(@times, abs(cellfun(@(x) x(2).qOs, wellSols)), dt)));
hold on;
plot(time,oip(1)-oip,'-.','LineWidth',3);
plot([0 1.2],oip([1 1]),'-k',time,min(t*oip(1),oip(1)),'-k', ...
    time(W(2).bt+[0 0]),[0 oip(1)],'--k');
hold off; axis tight; axis([0 max(time) 0 1.05*oip(1)]);
 xlabel('Time [days]','FontSize',14) 
end
ylabel('Cumulative Oil Production [m^3/day]','FontSize',14)


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
              
            end
    end
         
end
    

end
end







