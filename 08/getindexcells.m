clear Indx i j k nl
numc = 0;
ncel = deal(ny/8);
nl = 0;
nlay =3;
for i = 1: nlay: nx
    
   
for j = 1 : G.cells.num
     
         for k=1:nlay
   if G.cells.centroids(j,2) == nl + k - 0.5
         numc = numc + 1;
        Indx(numc) = G.cells.indexMap(j);
    end
    
end
end
 nl =  nl+k+nlay;
end
