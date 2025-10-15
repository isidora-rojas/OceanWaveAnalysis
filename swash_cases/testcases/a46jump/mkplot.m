% plots for hydraulic jump and drop test
%
load a46jmp02.tbl
fignr=7;
x=a46jmp02(:,1);
d=-a46jmp02(:,2);
h=a46jmp02(:,3);
s=a46jmp02(:,4);
plot(x,s,'k','LineWidth',1.0);hold on
plot(x,d,'k','LineWidth',2.0)
xlabel('x [m]');
ylabel('water depth [m]');
axis([0. 30.5 -0.5 0.3])
set(gca,'PlotBoxAspectRatio',[2 1 1]);
%orient tall;
eval(['print -dpng jump_' num2str(fignr) '.png']);
