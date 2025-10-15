% plots for Boers 1C experiment (Hm0, Tm02)
%
load l31set03
fignr=3;
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
%------------------------
smoot=1;
hm0=[];tm=[];
for k=1:mxc
   clear data df P F dof out
   data=wlv(:,k);
   nlen=length(data);
   if mod(nlen,2) == 1
      nlen=nlen-1;
   end
   [P,F,dof]=crosgk(data(1:nlen),data(1:nlen),nlen,smoot,dt,1,0);
   out.E(:,1) = real(P(:,1));
   out.f= F;
   out.dof = dof;
   df=out.f(2)-out.f(1);
   m0=df*sum(out.E(:,1));
   m2=df*sum(out.f(:).*out.f(:).*out.E(:,1));
   hm0=[hm0 4*sqrt(m0)];
   tm=[tm sqrt(m0/m2)];
end;
hm0=hm0';tm=tm';
H=subplot(1,2,1);
P1=plot(x,hm0);hold on;
%title('significant wave height')
axis([0 30 0.04 0.14])
xlabel('x [m]')
ylabel('H_{m0} [m]');
set(gca,'PlotBoxAspectRatio',[2 1 1]);
H=subplot(1,2,2);
P1=plot(x,tm);hold on;
%title('mean zero-crossing period')
axis([0 30 0.5 2.5])
xlabel('x [m]')
ylabel('T_{m02} [s]');
set(gca,'PlotBoxAspectRatio',[2 1 1]);
eval(['print -dpng boers_' num2str(fignr) '.png']);
