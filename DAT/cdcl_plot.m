clear all
load cdcl.dat

ll=1;
ms=1;

t=cdcl(:,1);

mb=2*ll;
me=mb+1;
cxcy=cdcl(:,mb:me);

mb=2*ms+2*ll;
me=mb+1;
xcyc=cdcl(:,mb:me);

mb=4*ms+2*ll;
me=mb+1;
xctyct=cdcl(:,mb:me);

mb=6*ms+2*ll;
me=mb+1;
thet=cdcl(:,mb:me);

figure(1)
plot(t,cxcy(:,1),'r-',t,cxcy(:,2),'b-')
xlabel('t')
ylabel('cd,cl')

figure(2)
plot(t,xcyc(:,1),'r-',t,xcyc(:,2),'b-')
xlabel('t')
ylabel('xc,yc')

figure(3)
plot(t,xctyct(:,1),'r-',t,xctyct(:,2),'b-')
xlabel('t')
ylabel('xct,yct')

figure(4)
plot(t,thet(:,1),'r-',t,thet(:,2),'b-')
xlabel('t')
ylabel('theta,thetat')


