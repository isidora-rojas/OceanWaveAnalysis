% plots for dam break test with wet bed everywhere
%
load a61dam02.tbl
fignr=2;
x=a61dam02(:,1);
wlv=a61dam02(:,2);
u=a61dam02(:,3);
subplot(2,1,1)
plot(x,wlv,'k');hold on
xlabel('x [m]');
ylabel('elevation [m]');
axis([0. 100. 0 1.1])
set(gca,'PlotBoxAspectRatio',[2 1 1]);
subplot(2,1,2)
plot(x,u,'k');hold on
xlabel('x [m]');
ylabel(' velocity [m/s]');
axis([0. 100. 0 3])
set(gca,'PlotBoxAspectRatio',[2 1 1]);
eval(['print -dpng dambr_' num2str(fignr) '.png']);
