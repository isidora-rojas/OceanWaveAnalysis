% plots for wind setup in semi-closed basin
%
load a34wnd01.tbl
t=a34wnd01(:,1)/3600;
wl=a34wnd01(:,2);
plot(t,wl,'k');hold;
load a34wnd02.tbl
t=a34wnd02(:,1)/3600;
wl=a34wnd02(:,2);
plot(t,wl,'k--');
load a34wnd03.tbl
t=a34wnd03(:,1)/3600;
wl=a34wnd03(:,2);
plot(t,wl,'k-.');
legend('\alpha=1','\alpha=0.5','\alpha=0.2')
xlabel('time [hr]');
ylabel('water level [m]');
set(gca,'PlotBoxAspectRatio',[2 1 1]);
print -dpng wind4.png
