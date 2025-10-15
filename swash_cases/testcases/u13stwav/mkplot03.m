% plots of vertical distribution of velocities for closed basin
%
fignr=10;
filname='u13stw02.tab';
kmax=10;                % number of layers
np=1;                   % number of locations
dt = 0.01;              % time step
t = [2.7; 3.6; 4.5];    % selected times
jz = 2;                 % layer interfaces
ju = 3;                 % U-velocity
jw = 5;                 % w-velocity
ml = ['r-';'g-';'b-'];  % markerline for model results
mle= ['ro';'go';'bo'];  % markerline for analytical solution
%--------------------------------------
load exactu1.dat
load exactu2.dat
load exactu3.dat
zu(1,:)=exactu1(:,1); ue(1,:)=exactu1(:,2);
zu(2,:)=exactu1(:,1); ue(2,:)=exactu2(:,2);
zu(3,:)=exactu1(:,1); ue(3,:)=exactu3(:,2);

load exactw1.dat
load exactw2.dat
load exactw3.dat
zw(1,:)=exactw1(:,1); we(1,:)=exactw1(:,2);
zw(2,:)=exactw1(:,1); we(2,:)=exactw2(:,2);
zw(3,:)=exactw1(:,1); we(3,:)=exactw3(:,2);

TMP = swash_loadTableVData(filname,np);

names = TMP{1};
unit  = TMP{2};
times = TMP{3};

for i=1:length(t)
    it(i) = find(times<t(i)+dt & times>t(i)-dt);
end

subplot(1,2,1)
for i=1:length(it)
    for k=1:kmax
        z(k) = 0.5*(TMP{2+jz,k}(it(i),1)+TMP{2+jz,k+1}(it(i),1));
        u(k) = TMP{2+ju,k}(it(i),1);
    end
    plot(u,z,ml(i,:)); hold on
end
xlabel(strcat('u-velocity [',strtrim(unit(ju)),']'));
ylabel(strcat('z [',strtrim(unit(jz)),']'));
axis([-0.2 0.2 -10 0])
legend(strcat('t=',num2str(t(1)),' s'),strcat('t=',num2str(t(2)),' s'),strcat('t=',num2str(t(3)),' s'),'Location','SouthEast')
for i=1:length(it)
    plot(ue(i,:),zu(i,:),mle(i,:));
end

subplot(1,2,2)
for i=1:length(it)
    for k=1:kmax+1
        z(k) = TMP{2+jz,k}(it(i),1);
        w(k) = TMP{2+jw,k}(it(i),1);
    end
    plot(w,z,ml(i,:)); hold on
    plot(we(i,:),zw(i,:),mle(i,:));
end
xlabel(strcat('w-velocity [',strtrim(unit(jw)),']'));
ylabel(strcat('z [',strtrim(unit(jz)),']'));
axis([-0.2 0.2 -10 0.1])

%orient tall;
suptitle(strcat('kh=\pi; 10 equidistant layers (standard)'));
eval(['print -dpng seiche' num2str(fignr) '.png']);
