    load shapet.dat
    load pt.dat
    load uct.dat
    load vct.dat
    load wot.dat
    load pht.dat

    [n,m]=size(pt);
    kn=n/nc;

    icontour=1;
    iquiver=0;
    label=0;
    for k=1:kn
      ks=(k-1)*nc+1;
      ke=ks+nc-1;
      if icontour==1
        cs=contour(xc,yc,wot(ks:ke,:),level);
        hold on
      end  
      if iquiver==1
        quiver(xc,yc,uct(ks:ke,:),vct(ks:ke,:))
        hold on
      end
      if label==1
        clabel(cs)
      end
      plot(shape0(:,2),shape0(:,3),'r-',shapee(:,2),shapee(:,3),'g-')
      hold on    
      ks=(k-1)*ns+1;
      ke=ks+ns-1;
      plot(st(ks:ke,3),st(ks:ke,4),'k-')
      hold off
      axis equal
      M(k)=getframe(gcf);
    end

    movie(M,5)
