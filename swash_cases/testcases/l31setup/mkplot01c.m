% spectra plots for Boers 1C experiment
%
load l31set03
fignr=1;
mxc=65;
dx=0.5;
dt=0.05;
%--------------------------------------
vars=whos;
j=0;
for i=1:length(vars)
    if ( strcmp(vars(i).name(1),'W') ),
       eval(['zeta = ' vars(i).name ';']);
       timestr = vars(i).name(8:17);
       [y,mo,d,h,m,s,ms]=str2dat(timestr,'hhttss_msc');
       j=j+1;
       time(j)=h*3600+m*60+s+ms/1000;
       wlv(j,:)=zeta;
    end;
end;
time=time';
%
x=Xp';
d=-Botlev';
%
aantal=8;
dis=[5 10 16 20 22 24 26 28];
emax=[0.005 0.005 0.004 0.004 0.003 0.003 0.002 0.001];
smoot=100;
flow=0;
fhigh=2.0;
%------------------------
sub=1;
for i=1:aantal
   clear data P F dof out
   j=dis(i)/dx+1;
   data=wlv(:,j);
   nlen=length(data);
   if mod(nlen,2) == 1
      nlen=nlen-1;
   end
   [P,F,dof]=crosgk(data(1:nlen),data(1:nlen),nlen,smoot,dt,1,0);
   out.E(:,1) = real(P(:,1));
   out.f= F;
   out.dof = dof;
   H=subplot(aantal/2,2,sub);
   P1=plot(out.f(:),out.E(:,1));hold on;
   clear data P F dof out
   title(['x = ' num2str(dis(i)) ' m'])
   if i==aantal-1 | i==aantal
      xlabel('f [Hz]')
   end
   if mod(i,2)==1
      ylabel('E [m^2/Hz]')
   end
   axis([flow fhigh 0 emax(i)])
   set(gca,'PlotBoxAspectRatio',[2 1 1]);
   sub=sub+1;
end
orient tall;
eval(['print -dpng boers_' num2str(fignr) '.png']);
