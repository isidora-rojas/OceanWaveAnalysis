% plots for wind setup in closed basin
%
load a31wnd01.tbl
t=a31wnd01(:,1)/3600;
wl=a31wnd01(:,2);
plot(t,wl,'k');hold;
load a31wnd02.tbl
t=a31wnd02(:,1)/3600;
wl=a31wnd02(:,2);
plot(t,wl,'k--');
legend('absolute','relative')
xlabel('time [hr]');
ylabel('water level [m]');
set(gca,'PlotBoxAspectRatio',[2 1 1]);
print -dpng wind1.png
