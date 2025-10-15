% spectra plotjes tbv koppeling SWAN - SWASH
%
clear wlv
fignr=2;
nloc=1;
dt=0.1;
%--------------------------------------
TMP = swash_loadTableData('t11cou02.tbl',nloc,6);
time = TMP{2};
x    = TMP{3}(1,:);
y    = TMP{4}(1,:);
wlv  = TMP{7};
wlv(wlv==-99)=NaN;
%
smoot=30;
%------------------------
load t11cou02.sp1
swans=t11cou02;
swans(swans==-99)=NaN;
%
clear data P F dof out
data=wlv(:,1);
nlen=length(data);
if mod(nlen,2) == 1
   nlen=nlen-1;
end
[P,F,dof]=crosgk(data(1:nlen),data(1:nlen),nlen,smoot,dt,2,0);
out.E(:,1) = real(P(:,1));
out.f= F;
out.dof = dof;
P1=plot(out.f(:),out.E(:,1),'k');hold on;
P2=plot(swans(:,1),swans(:,2),'k');
set(P2,'linewidth',1.3)
axis([0 0.5 0 0.14])
legend('SWASH','SWAN')
xlabel('f [Hz]')
ylabel('E [m^2/Hz]')
eval(['print -dpng coupling' num2str(fignr) '.png']);
