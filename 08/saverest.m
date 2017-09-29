function saverest(dir1,file,def,pod,dpod,ts,dv,preport,training,pri)



text = [dir1 file];
for i=1:ts
    its(i,1)=preport(1,i).iterations;
end

if ~exist(text, 'file')
    fileID = fopen(text,'w');
    fprintf(fileID,'\\begin{table}[!ht]\\centering\n');
    fprintf(fileID,'\\begin{minipage}{0.5\\textwidth}\n');
    fprintf(fileID,' \\centering\n');
    fprintf(fileID,'\\begin{tabular}{ ||c|c|c|c||} \n');
    fprintf(fileID,'\\hline\n');
    fprintf(fileID,'Total&Method &Iterations &\\%% of ICCG\\\\ \n');
    fprintf(fileID,'ICCG& & &Iterations\\\\ \n');
    fclose(fileID);
        
end

if (def==0)&&(training==1)
    ttits=sum(its)
    %         fileID = fopen(text,'a');
    %         fprintf(fileID,'\\hline  \n');
    %         fprintf(fileID,['$10^{' num2str(per) '}$ &' num2str(ttits) '& DICCG$_{' num2str(dv) '}$& & &  & \\\\ \n']);
    %         fclose(fileID);
    save([dir1 num2str(pri)  'ttits.mat'],'ttits')
    
else if (def==1)&&(training==1)
    tits=sum(its);
%    icits=sum(its(1:dv));
%    dicits=sum(its(dv+1:ts));
    load([dir1  num2str(pri)  'ttits.mat'],'ttits')
    perc=round(mean(tits*100/ttits));
    if pod==0
        
        fileID = fopen(text,'a');
        fprintf(fileID,'\\hline \n');
        fprintf(fileID,[ num2str(ttits) '& DICCG$_{' num2str(dv) '}$&' num2str(tits)  '&' num2str(perc) '\\\\ \n']);
        fclose(fileID);
    else
        fileID = fopen(text,'a');
        fprintf(fileID,'\\hline  \n');
        fprintf(fileID,[ num2str(ttits) '& DICCG$_{POD_{' num2str(numel(dpod)) '}}$&' num2str(tits)  '&' num2str(perc) ' \\\\ \n']);
        fclose(fileID);
        if numel(dpod)==10
            fileID = fopen(text,'a');
            fprintf(fileID,'\\hline  \n');
            fprintf(fileID,'\\end{tabular} \n');
            fprintf(fileID,'\\caption{Comparison between the ICCC and DICCG methods of the average number of linear iterations. }\\label{table:litertot} \n');
            fprintf(fileID,'\\end{minipage}  \n');
            fprintf(fileID,'\\end{table}  \n');
            fclose(fileID);
        end
    end
    
    end
end
end

