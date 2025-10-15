% plots for short wave in closed basin
%
load exact.dat
load a13stw03.tbl
t=a13stw03(:,1);
wl=a13stw03(:,2);
plot(t,wl,'k');hold;
plot(exact(1:4:601,1),exact(1:4:601,2),'ko');
xlabel('time [s]');
ylabel('water level [m]');
title('kh=10\pi; 3 layers: 5%, 25% 70%')
set(gca,'PlotBoxAspectRatio',[2 1 1]);
print -dpng seiche7.png
