function addplots(dir1,dir2,file,file1)
for i = 1 : numel(file)
text1 = [dir1 file1];
text = [dir2 file{i}];
if ~exist(text1, 'file')
    fileID = fopen(text,'w');
    fprintf(fileID,'\\begin{figure}[!h] \\hspace{-1cm}\n');
    fprintf(fileID,'\\begin{minipage}{0.9\\textwidth}\n');
    fprintf(fileID,['\\includegraphics[width=4cm,height=4cm,keepaspectratio]{' text '.jpg} \n']);  
    fprintf(fileID,['\\caption{' file{i} '}\n']);
    fprintf(fileID,['\\label{fig:' file{i} '}\n']);
    fprintf(fileID,'\\end{minipage}%% \n');
   fprintf(fileID,'\\end{figure} \n');
    fclose(fileID);
else
        fileID = fopen(text,'a');
    fprintf(fileID,'\\begin{figure}[!h] \\hspace{-1cm}\n');
    fprintf(fileID,'\\begin{minipage}{0.9\\textwidth}\n');
    fprintf(fileID,['\\includegraphics[width=4cm,height=4cm,keepaspectratio]{' text '.jpg} \n']);  
    fprintf(fileID,['\\caption{' file{i} '}\n']);
    fprintf(fileID,['\\label{fig:' file{i} '}\n']);
    fprintf(fileID,'\\end{minipage}%% \n');
   fprintf(fileID,'\\end{figure} \n');
    fclose(fileID);

end
end
end