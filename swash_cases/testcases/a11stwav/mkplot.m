% plots for short wave in closed basin
%
load exact.dat
load a11stw01.tbl
t=a11stw01(:,1);
wl=a11stw01(:,2);
plot(t,wl,'k');hold;
plot(exact(1:4:601,1),exact(1:4:601,2),'ko');
xlabel('time [s]');
ylabel('water level [m]');
title('kh=0.55\pi; depth averaged')
set(gca,'PlotBoxAspectRatio',[2 1 1]);
print -dpng seiche1.png
