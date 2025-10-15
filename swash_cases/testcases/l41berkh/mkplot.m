% plots for wave over elliptic shoal
%
clear all
fignr=1;
mxc=381;
myc=221;
%--------------------------------------
TMP = swash_loadTableData('section1.tbl',mxc,5);
%t   = TMP{2};
x   = TMP{3}(1,:);x=x';
wlv = TMP{4};
wlv(wlv==-99)=NaN;
%
for i=1:mxc
    h(i)=max(wlv(:,i))-min(wlv(:,i));
end
h=h';h=h/.0464;
subplot(4,2,1)
plot(x,h);hold on
title('section 1')
axis([-5 5 0 2.5])
xlabel('x [m]')
ylabel('H/H_0')
clear TMP wlv h
TMP = swash_loadTableData('section2.tbl',mxc,5);
wlv = TMP{4};
wlv(wlv==-99)=NaN;
%
for i=1:mxc
    h(i)=max(wlv(:,i))-min(wlv(:,i));
end
h=h';h=h/.0464;
subplot(4,2,2)
plot(x,h);hold on
title('section 2')
axis([-5 5 0 2.5])
xlabel('x [m]')
clear TMP wlv h
TMP = swash_loadTableData('section3.tbl',mxc,5);
wlv = TMP{4};
wlv(wlv==-99)=NaN;
%
for i=1:mxc
    h(i)=max(wlv(:,i))-min(wlv(:,i));
end
h=h';h=h/.0464;
subplot(4,2,3)
plot(x,h);hold on
title('section 3')
axis([-5 5 0 2.5])
xlabel('x [m]')
ylabel('H/H_0')
clear TMP wlv h
TMP = swash_loadTableData('section4.tbl',mxc,5);
wlv = TMP{4};
wlv(wlv==-99)=NaN;
%
for i=1:mxc
    h(i)=max(wlv(:,i))-min(wlv(:,i));
end
h=h';h=h/.0464;
subplot(4,2,4)
plot(x,h);hold on
title('section 4')
axis([-5 5 0 2.5])
xlabel('x [m]')
clear TMP wlv h
TMP = swash_loadTableData('section5.tbl',mxc,5);
wlv = TMP{4};
wlv(wlv==-99)=NaN;
%
for i=1:mxc
    h(i)=max(wlv(:,i))-min(wlv(:,i));
end
h=h';h=h/.0464;
subplot(4,2,5)
plot(x,h);hold on
title('section 5')
axis([-5 5 0 2.5])
xlabel('x [m]')
ylabel('H/H_0')
clear TMP x wlv h
TMP = swash_loadTableData('section6.tbl',myc,5);
x   = TMP{3}(1,:);x=x';
wlv = TMP{4};
wlv(wlv==-99)=NaN;
%
for i=1:myc
    h(i)=max(wlv(:,i))-min(wlv(:,i));
end
h=h';h=h/.0464;
subplot(4,2,6)
plot(x,h);hold on
title('section 6')
axis([0 10 0 2.5])
xlabel('y [m]')
clear TMP wlv h
TMP = swash_loadTableData('section7.tbl',myc,5);
wlv = TMP{4};
wlv(wlv==-99)=NaN;
%
for i=1:myc
    h(i)=max(wlv(:,i))-min(wlv(:,i));
end
h=h';h=h/.0464;
subplot(4,2,7)
plot(x,h);hold on
title('section 7')
axis([0 10 0 2.5])
xlabel('y [m]')
ylabel('H/H_0')
clear TMP wlv h
TMP = swash_loadTableData('section8.tbl',myc,5);
wlv = TMP{4};
wlv(wlv==-99)=NaN;
%
for i=1:myc
    h(i)=max(wlv(:,i))-min(wlv(:,i));
end
h=h';h=h/.0464;
subplot(4,2,8)
plot(x,h);hold on
title('section 8')
axis([0 10 0 2.5])
xlabel('y [m]')
clear TMP x y wlv h
orient tall;
eval(['print -dpng berkhoff_' num2str(fignr) '.png']);
