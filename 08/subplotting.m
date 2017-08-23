clim=[-1 2]; % Data range...

figure; hold on; axis off;
ax(1)=axes('position',[0.05 0.55 0.35 0.40]);
  imagesc(rand(5,5)-1);
ax(2)=axes('position',[0.05 0.05 0.35 0.40]);
  imagesc(rand(6,6));
ax(3)=axes('position',[0.45 0.55 0.35 0.40]);
  imagesc(rand(7,7));
ax(4)=axes('position',[0.45 0.05 0.35 0.40]);
  imagesc(2*rand(8,8));

caxis(clim);
set(ax,'CLim',clim);
h=colorbar;
  set(h,'position',[0.85 0.05 0.10 0.90])