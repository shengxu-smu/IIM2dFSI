clear all

load surface1.dat

theta=180/3.1415926;

ll=1;

if ll==1
    f=surface1;
elseif ll==2
    f=surface2;
elseif ll==3
    f=surface3;
elseif ll==4
    f=surface4;
elseif ll==5
    f=surface5;
end

figure(1)
plot(f(:,2),f(:,3),'ko',f(:,4),f(:,5),'g-')
xlabel('xs')
ylabel('ys')
axis equal

figure(2)
plot(theta*f(:,1),f(:,10),'ko',theta*f(:,1),f(:,12),'g-')
xlabel('alpha')
ylabel('us')

figure(3)
plot(theta*f(:,1),f(:,11),'ko',theta*f(:,1),f(:,13),'g-')
xlabel('alpha')
ylabel('vs')

figure(4)
plot(theta*f(:,1),f(:,14),'-o',theta*f(:,1),f(:,15),'--')
xlabel('alpha')
ylabel('fn,ft')
 
ip1=3;
ip2=5;
ip=20+(ip1-1)+(ip2-1)*3;

figure(5)
plot(theta*f(:,1),f(:,ip),'-')
xlabel('alpha')
ylabel('jc')
