% plots for linear progressive wave in a flume (kh = 12.57)
%
a0=0.04;
T=3.5791;
tt=0:0.01:72;
w=2*pi/T;
y=cos(w*tt);
load a14prw04.tbl
t=a14prw04(:,1)/T;
wl=a14prw04(:,2)/a0;
plot(t(1:10:5760),wl(1:10:5760),'ko','MarkerSize',2);hold;
plot(tt/T,y,'k');
axis([10 20 -2 2])
xlabel('t/T_0');
ylabel('\zeta/a_0');
title('kh=4\pi; 3 layers: 10%, 20% 70%')
set(gca,'PlotBoxAspectRatio',[2 1 1]);
print -dpng prowave4.png
