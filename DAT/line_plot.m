load gx.dat;
load gy.dat;
ip=2:5;

figure(1)
plot(gx(:,1),gx(:,ip),'--')
legend('u','v','p','o')

figure(2)
plot(gy(:,1),gy(:,ip),'--')
legend('u','v','p','o')

load cdcl.dat;
[m,n]=size(cdcl);
k=1;

figure(3)
plot(cdcl(1:k:m,1),cdcl(1:k:m,2),'r--')
xlabel('t')
ylabel('cd')

figure(4)
plot(cdcl(1:k:m,1),cdcl(1:k:m,3),'r--')
xlabel('t')
ylabel('cl')

load rangle.dat;

figure(5)
polar(rangle(:,1),rangle(:,2),'-')
