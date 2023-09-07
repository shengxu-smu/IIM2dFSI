clear all

load xc.dat
load yc.dat
load xe.dat
load ye.dat

load u.dat
load v.dat
load p.dat
load d.dat

load eu.dat
load ev.dat
load ep.dat

mc=size(xc,1);
nc=size(yc,1);
me=size(xe,1);
ne=size(ye,1);

ixes=1;
ixee=me;
ixcs=1;
ixce=mc;
jyes=1;
jyee=ne;
jycs=1;
jyce=nc;

figure(1)
surf(xe(ixes:ixee),yc(jycs:jyce),u(jycs:jyce,ixes:ixee))
figure(2)
surf(xe(ixes:ixee),yc(jycs:jyce),eu(jycs:jyce,ixes:ixee))

figure(3)
surf(xc(ixcs:ixce),ye(jyes:jyee),v(jyes:jyee,ixcs:ixce))
figure(4)
surf(xc(ixcs:ixce),ye(jyes:jyee),ev(jyes:jyee,ixcs:ixce))

figure(5)
surf(xc(ixcs:ixce),yc(jycs:jyce),p(jycs:jyce,ixcs:ixce))
figure(6)
surf(xc(ixcs:ixce),yc(jycs:jyce),ep(jycs:jyce,ixcs:ixce))

figure(7)
surf(xc(ixcs:ixce),yc(jycs:jyce),d(jycs:jyce,ixcs:ixce))

