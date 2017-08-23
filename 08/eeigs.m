for cp = [0 1]
    
    for per=[3:2:7]
        
        for def=[1 ]
            
            % Number of deflation vectors
            
            % If we want to use POD pod==1
            for pod = [0 1]
                clearvars -except per def pod cp
                
                
                dpod=[];
                if (def == 0)  && (pod == 1)
                    pod=0;
                    continue
                end
                for rr=1:2
                    if (pod == 0) && (rr==2)
                        continue
                    end
                    vars
                    close all
                    dpod=podv{rr};
                    if pod== 0
                        dpod=[];
                    end
                    
                    
                    
                    %Create the directory
                    dir='/mnt/sda2/cortes/Results/17_05/two_phases/eigs/08/';
                    
                    folder=[ '10-' num2str(k) '_' num2str(sz) 'nz' num2str(nz) 'perm_' num2str(per) 'cp' num2str(cp)];
                    %mkdir([dir], folder)
                    dir1 = [dir folder '/'];
                    
                    folder=[  'def_' num2str(def) '_pod_' num2str(numel(dpod)) ];
                    %mkdir([dir1], folder)
                    dir2 = [dir1 folder '/'];
                        B=[dir2    'preport.mat'];
                    load(B)
                    
                    ts = length(preport);
                    for i = 11 : ts
                        cnn(i) = preport(1,i).cond(2);
                    end
                    
                    
                    f(11) = figure(11);
                    
                    hn=plot((11*dT:dT:T)/day,cnn(11:ts),'b*');
                    ylim([0 1e29])
                    % legend('ICCG');
                    axis square
                    ylabel(' \kappa_2(E)','FontSize',16)
                    xlabel('Time (days) ','FontSize',16)
                    file{11} = 'econdest_s';
                    B=[dir2   file{11}   '.jpg'];
                    
                       saveas(f(11),B)
                        B=[dir2   file{11}   '.m'];
                    
                       saveas(f(11),B)
                    %           f(10) = figure(10);
                    %             hn=plot((dT:dT:10*dT)/day,cnn(1:10),'r*');
                    %             hold on
                    %             hn=plot((11*dT:dT:T)/day,cnn(11:ts),'b*');
                    %        % legend('ICCG');
                    %         axis square
                    %     ylabel(' \kappa_2','FontSize',16)
                    %     xlabel('Time (days) ','FontSize',16)
                    %     file{10} = 'condest_1';
                    %    B=[dir2   file{10}   '.jpg'];
                    %  saveas(f(10),B)
                    
                    
                    
                end
            end
        end
    end
end