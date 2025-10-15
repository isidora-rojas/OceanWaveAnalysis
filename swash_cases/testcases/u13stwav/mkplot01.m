% plots for short wave in closed basin
%
load exact.dat
load u13stw01.tbl
t=u13stw01(:,1);
wl=u13stw01(:,2);
plot(t,wl,'k');hold;
plot(exact(1:4:601,1),exact(1:4:601,2),'ko');
xlabel('time [s]');
ylabel('water level [m]');
title('kh=\pi; 2 equidistant layers (box)')
set(gca,'PlotBoxAspectRatio',[2 1 1]);
print -dpng seiche8.png
