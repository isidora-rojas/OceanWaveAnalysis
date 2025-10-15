% plots for tidal wave flow
%
load a81tid01.tbl
fignr=1;
x=a81tid01(:,1);
d=a81tid01(:,2);
s=a81tid01(:,3);
u=a81tid01(:,4);
subplot(2,1,1)
plot(x,s,'k','LineWidth',1.5);hold on
plot(x,-d,'k')
z=x-x+2.18;
plot(x,z,'ro')
xlabel('x [m]');
ylabel('elevation [m]');
axis([0. 14000. -60 10])
set(gca,'PlotBoxAspectRatio',[2 1 1]);
subplot(2,1,2)
plot(x,u,'k','LineWidth',1.5);hold on
load exact.u
plot(exact(1:28:1400,1),exact(1:28:1400,2),'ro')
xlabel('x [m]');
ylabel(' velocity [m/s]');
axis([0. 14000. 0 0.2])
set(gca,'PlotBoxAspectRatio',[2 1 1]);
eval(['print -dpng tidal_' num2str(fignr) '.png']);
