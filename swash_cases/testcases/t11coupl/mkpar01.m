% results coupling SWAN - SWASH (Hm0, Tm01, Tm02)
%
clear wlv
dt=0.1;
nloc=1;
is=121;
ie=1501;
%
%--------------------------------------
TMP = swash_loadTableData('t11cou01.tbl',nloc,6);
time = TMP{2};
x    = TMP{3}(1,:);
y    = TMP{4}(1,:);
wlv  = TMP{7};
wlv(wlv==-99)=NaN;
%
load t11cou01.tab
hm0s =t11cou01(:,5);
tm01s=t11cou01(:,7);
tm02s=t11cou01(:,8);
%
%------------------------
smoot=1;
hm0=[];tm01=[];tm02=[];
for k=1:nloc
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
   out.f(is)
   out.f(ie)
   m0=df*sum(out.E(is:ie,1));
   m1=df*sum(out.f(is:ie).*out.E(is:ie,1));
   m2=df*sum(out.f(is:ie).*out.f(is:ie).*out.E(is:ie,1));
   hm0=[hm0 4*sqrt(m0)];
   tm01=[tm01 m0/m1];
   tm02=[tm02 sqrt(m0/m2)];
end;
%
h=[hm0s hm0; tm01s tm01; tm02s tm02]
