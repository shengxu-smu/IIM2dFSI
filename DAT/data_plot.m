clear all

load xc.dat
load yc.dat

load uc.dat
load vc.dat
load p.dat
load d.dat
load wo.dat
load ph.dat

load surface1.dat

kip=2;

icontour=1;
ishape=1;
iquiver=1;

igiven=0;
label=0;
level=60;

mc=size(xc,1);
nc=size(yc,1);

ixcs=1;
ixce=mc;
jycs=1;
jyce=nc;

figure(1)
surf(xc(ixcs:kip:ixce),yc(jycs:kip:jyce),uc(jycs:kip:jyce,ixcs:kip:ixce))

figure(2)
surf(xc(ixcs:kip:ixce),yc(jycs:kip:jyce),vc(jycs:kip:jyce,ixcs:kip:ixce))

figure(3)
surf(xc(ixcs:kip:ixce),yc(jycs:kip:jyce),p(jycs:kip:jyce,ixcs:kip:ixce))

figure(4)
surf(xc(ixcs:kip:ixce),yc(jycs:kip:jyce),d(jycs:kip:jyce,ixcs:kip:ixce))

wov=[-2.5:0.1:2.5];
phv1=[-4:0.2:4];
phv2=[-0.1:0.01:0.1];
phv=[phv1,phv2];
pv=[-4:0.2:4];

iplot=1;

figure(5)
if icontour==1
  if iplot==1
    if igiven==0
      H=pcolor(xc,yc,wo);
      shading interp;
      caxis([-3 3]);
      axis equal;
    else
      cs=contour(xc,yc,wo,wov);
    end
  end
  if iplot==2
    if igiven==0
      cs=contour(xc,yc,ph,level);
    else
      cs=contour(xc,yc,ph,phv);
    end
  end
  if iplot==3
    if igiven==0
      cs=contour(xc,yc,p,level);
    else
      cs=contour(xc,yc,p,pv);
    end
  end  
  if label==1
    clabel(cs)
  end
  hold on
end  
if iquiver==1
  quiver(xc(ixcs:kip:ixce),yc(jycs:kip:jyce),...
         uc(jycs:kip:jyce,ixcs:kip:ixce),vc(jycs:kip:jyce,ixcs:kip:ixce))
  hold on
end
if ishape==1
  plot(surface1(:,2),surface1(:,3),'k-')
end
hold off
axis equal

