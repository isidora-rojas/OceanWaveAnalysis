% plots for linear progressive wave in a flume (kh = 3.14)
%
a0=0.01;
T=3.5791;
tt=0:0.01:72;
w=2*pi/T;
y=cos(w*tt);
load a14prw01.tbl
t=a14prw01(:,1)/T;
wl=a14prw01(:,2)/a0;
plot(t(1:5:2881),wl(1:5:2881),'ko','MarkerSize',2);hold;
plot(tt/T,y,'k');
axis([10 20 -2 2])
xlabel('t/T_0');
ylabel('\zeta/a_0');
title('kh=\pi; 2 equidistant layers')
set(gca,'PlotBoxAspectRatio',[2 1 1]);
print -dpng prowave1.png
