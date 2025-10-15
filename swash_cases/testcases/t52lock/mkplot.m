% plots of vertical distribution of salinity for lock exchange test
%
kmax=30;
%--------------------------------------
%
clear x zk
load t52lok01
%
np=length(Xp);
%
for j=1:np
    for k=1:kmax
        if k < 10
           eval(['s=Salk0' num2str(k) '_000015_000;'])
           eval(['z=zk0' num2str(k) '_000015_000;'])
        else
           eval(['s=Salk' num2str(k) '_000015_000;'])
           eval(['z=zk' num2str(k) '_000015_000;'])
        end
        sal(j,k)=s(j);
        zk(j,k)=z(j);
        x(j,k)=Xp(j);
    end
end
subplot(2,1,1)
pcolor(x,zk+0.3,double(sal));shading flat;axis equal
axis([0.01 2 0.0 0.29])
caxis([-0.5 5.5]);
title('hydrostatic mode')
ylabel('z (m)');
%
clear s sal
load t52lok02
%
for j=1:np
    for k=1:kmax
        if k < 10
           eval(['s=Salk0' num2str(k) '_000015_000;'])
           eval(['z=zk0' num2str(k) '_000015_000;'])
        else
           eval(['s=Salk' num2str(k) '_000015_000;'])
           eval(['z=zk' num2str(k) '_000015_000;'])
        end
        sal(j,k)=s(j);
        zk(j,k)=z(j);
        x(j,k)=Xp(j);
    end
end
subplot(2,1,2)
pcolor(x,zk+0.3,double(sal));shading flat;axis equal
axis([0.01 2 0.0 0.29])
caxis([-0.5 5.5]);
title('non-hydrostatic mode')
ylabel('z (m)');
xlabel('x (m)');
print -dpng lock01.png
