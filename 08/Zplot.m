 load([dir1 'Z.mat'],'Z')  
                        size(Z)
figure
for i=1:30
    plot(Z(:,i))
     pause
end