%%If we want to save the files/graphs, it's neccesary to run this part
%                 
%                 file{1} = 'Permeability';
%                 file{2} = 'RelPerm';
%                 file{3} = 'Solution';
%                 file{4} = 'Iterations';
%                 file{5} = 'Est_Cond_number';
%                 file{6} = 'Cond_number';
%                 if cp==1
%                     
%                     filecp='cp';
%                     axis equal tight
%                     B=[dir2   filecp  '.fig'];
%                     saveas(hcp,B)
%                     B=[dir2   filecp   '.jpg'];
%                     saveas(hcp,B)
%                 end
%                 for i=1:numel(f)
%                    
%                     axis equal tight
%                     axis square
%                     %i
%                 B=[dir2   file{i}  '.fig'];
%                 savefig(f(i),B)
%                  B=[dir2   file{i}   '.jpg'];
%                  saveas(f(i),B)
%                 end
%filetx = ['results.txt'];
%if solver == [1 2 3 4]

               % saveres(dir1,filetx,def,pod,dpod,per,ts,dv,preport)
%end
                addplots(dir1,dir2,file,'plots.txt')
                clear f
                filews=['workspace'];
                filename=[dir2 filews];
                save(filename)