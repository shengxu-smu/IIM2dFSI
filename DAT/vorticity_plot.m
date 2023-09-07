clear all

load xc.dat
load yc.dat
load wo.dat

load surface1.dat
load surface2.dat
load surface3.dat
load surface4.dat
load surface5.dat

igiven=0;
wov=[-2.5:0.1:2.5];

figure(1)
if igiven==0
  H=pcolor(xc,yc,wo);
  shading interp;
  caxis([-4 4]);
  axis equal;
else
  cs=contour(xc,yc,wo,wov);
end
hold on
plot(surface1(:,2),surface1(:,3),'k-')
plot(surface2(:,2),surface2(:,3),'k-')
plot(surface3(:,2),surface3(:,3),'k-')
plot(surface4(:,2),surface4(:,3),'k-')
plot(surface5(:,2),surface5(:,3),'k-')
hold off
axis equal

