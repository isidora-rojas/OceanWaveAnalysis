% plots for Boers 1C experiment (Hm0, setup)
%
load l31set03.tab
fignr=2;
%------------------------
H=subplot(1,2,1);
P1=plot(l31set03(:,1),l31set03(:,2),'k');hold on;
title('significant wave height')
axis([0 30 0.04 0.14])
xlabel('x [m]')
ylabel('H_{m0} [m]');
set(gca,'PlotBoxAspectRatio',[2 1 1]);
H=subplot(1,2,2);
P1=plot(l31set03(:,1),1000*l31set03(:,3),'k');hold on;
title('wave-induced setup')
axis([0 30 -3 8])
xlabel('x [m]')
ylabel('E[\eta] [mm]');
set(gca,'PlotBoxAspectRatio',[2 1 1]);
eval(['print -dpng boers_' num2str(fignr) '.png']);
