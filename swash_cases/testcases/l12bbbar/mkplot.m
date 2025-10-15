% plots for regular wave over submerged bar
%
fignr=2;
%------------simulation----------------
load station4.tbl
t=station4(:,1);wlv=station4(:,2);
load station5.tbl
t=[t station5(:,1)];wlv=[wlv station5(:,2)];
load station6.tbl
t=[t station6(:,1)];wlv=[wlv station6(:,2)];
load station7.tbl
t=[t station7(:,1)];wlv=[wlv station7(:,2)];
load station8.tbl
t=[t station8(:,1)];wlv=[wlv station8(:,2)];
load station9.tbl
t=[t station9(:,1)];wlv=[wlv station9(:,2)];
load station10.tbl
t=[t station10(:,1)];wlv=[wlv station10(:,2)];
load station11.tbl
t=[t station11(:,1)];wlv=[wlv station11(:,2)];
%---------------------------------------
sti  = [' 4';' 5';' 6';' 7';' 8';' 9';'10';'11'];
dis  = ['10.5 m';'12.5 m';'13.5 m';'14.5 m';'15.7 m';'17.3 m';'19.0 m';'21.0 m'];
labl = ['105m';'125m';'135m';'145m';'157m';'173m';'190m';'210m'];
sub=1;
for i=1:8,
   H=subplot(4,2,sub);
   P1=plot(t(:,i)-0.15,100*wlv(:,i)); hold on;
   title(['station ' sti(i,:) ' (x = ' dis(i,:) ')'])
   axis([33 39 -2.0 4.0])
   if i==7 | i==8
      xlabel('time [s]');
   end
   if i==1 | i==3 | i==5 | i==7
      ylabel('water level [cm]');
   end
   set(gca,'PlotBoxAspectRatio',[2 1 1]);
   sub=sub+1;
end;
orient tall;
eval(['print -dpng bb_' num2str(fignr) '.png']);
