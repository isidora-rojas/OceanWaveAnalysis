% plots of vertical distribution of flow velocity and turbulence quantities
%
fignr=2;
casenr='u33tur02';
kmax=10;                % number of layers
np=1;                   % number of locations
jz = 2;                 % layer interfaces
ju = 3;                 % U-velocity
jk = 5;                 % turbulent kinetic energy
je = 6;                 % dissipation rate of turbulent energy
jv = 7;                 % vertical eddy viscosity
%--------------------------------------

d=0.1;
rho=1023;
mu=0.001;
nu=mu/rho;

close all
clear TMP z u

filname=[casenr '.tbl'];
eval(['load ' filname]);
eval(['h=' casenr '(1,1);']);
eval(['ustar=' casenr '(1,2);']);
ustar2=ustar*ustar;

filname=[casenr '.tab'];
TMP = swash_loadTableVData(filname,np);

subplot(2,2,1)
for k=1:kmax
    z(k) = 0.5*(TMP{2+jz,k}(1,1)+TMP{2+jz,k+1}(1,1));
    u(k) = TMP{2+ju,k}(1,1);
end
zp=(z+d)*ustar/nu;
up=u/ustar;
semilogx(zp,up,'k'); hold on
xlabel('z^+ [-]')
ylabel('u/u^* [-]')
axis([10 10000 10 30])

subplot(2,2,2)
for k=1:kmax+1
    z(k)   = TMP{2+jz,k}(1,1);
    tke(k) = TMP{2+jk,k}(1,1);
end
plot(tke/ustar2,(z+d)/h,'k'); hold on
xlabel('k/u^{*2} [-]')
ylabel('z/h [-]')
axis([0 5 0 1])

subplot(2,2,3)
for k=1:kmax+1
    z(k)   = TMP{2+jz,k}(1,1);
    eps(k) = TMP{2+je,k}(1,1);
end
plot(eps/(ustar*ustar2),(z+d)/h,'k'); hold on
xlabel('\epsilon/u^{*3} [-]')
ylabel('z/h [-]')
axis([0 2500 0 1])

subplot(2,2,4)
for k=1:kmax+1
    z(k)    = TMP{2+jz,k}(1,1);
    visc(k) = TMP{2+jv,k}(1,1);
end
plot(visc/(ustar*h),(z+d)/h,'k'); hold on
xlabel('\nu_t/hu^{*} [-]')
ylabel('z/h [-]')
axis([0 0.12 0 1])

eval(['print -dpng turbchannel' num2str(fignr) '.png']);
