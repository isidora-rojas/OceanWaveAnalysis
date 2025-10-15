% plots for linear progressive wave in a flume (kh = 15.71)
%
a0=0.05;
T=3.5791;
tt=0:0.01:72;
w=2*pi/T;
y=cos(w*tt);
load a14prw05.tbl
t=a14prw05(:,1)/T;
wl=a14prw05(:,2)/a0;
plot(t(1:10:5760),wl(1:10:5760),'ko','MarkerSize',2);hold;
plot(tt/T,y,'k');
axis([10 20 -2 2])
xlabel('t/T_0');
ylabel('\zeta/a_0');
title('kh=5\pi; 3 layers: 5%, 15% 80%')
set(gca,'PlotBoxAspectRatio',[2 1 1]);
print -dpng prowave5.png
