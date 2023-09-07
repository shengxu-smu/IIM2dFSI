c-----------------------------------------------------------------------
c
      subroutine jc_firstsecond
      include 'parameter.inc'
      include 'surface.inc'
      real*8 gacobi,gacobi2,dgacobi,dnx,dny,r1,r2,r3,ts,tm
      real*8 ft,fn,tx,ty,dft,dfn,dtx,dty,ddft,ddfn,gt,gn,g1,g2,dgt,dgn
      real*8 duxjc,duyjc,dduxjc,dduyjc,dduxyjc,ddpxjc0,ddpxjc1
      real*8 dvxjc,dvyjc,ddvxjc,ddvyjc,ddvxyjc,ddpyjc0,ddpyjc1
      real*8 dpxjc,dpyjc      

      DO l=1,ms

      do m=0,ns-1
        tx=taox(m,l)
        ty=taoy(m,l)
        dtx=dtaox(m,l)
        dty=dtaoy(m,l)
        ft=fx(m,l)
        fn=fy(m,l)
        dft=dfx(m,l)
        dfn=dfy(m,l)
        ddft=ddfx(m,l)
        ddfn=ddfy(m,l)

        gacobi=sqrt(tx*tx+ty*ty)
        gacobi2=tx*tx+ty*ty
        dgacobi=(tx*dtx+ty*dty)/gacobi
        dnx=dty/gacobi-ty*dgacobi/gacobi2
        dny=-dtx/gacobi+tx*dgacobi/gacobi2
        ts=(tx*tx-ty*ty)/gacobi2
        tm=2.0d0*tx*ty/gacobi

        gt=gx(m,l)
        gn=gy(m,l)
        g1=(gt*tx+gn*ty)/gacobi
        g2=(gt*ty-gn*tx)/gacobi
        dgt=dgx(m,l)
        dgn=dgy(m,l)

        duxjc=-Re*ty*tx*ft/gacobi2
        duyjc=Re*tx*tx*ft/gacobi2

        dvxjc=-Re*ty*ty*ft/gacobi2
        dvyjc=Re*tx*ty*ft/gacobi2

        dpxjc=(tx*dfn+ty*(dft+gacobi*gn))/gacobi2
        dpyjc=(ty*dfn-tx*(dft+gacobi*gn))/gacobi2

        r1=-dtx*duxjc-dty*duyjc
        r2=-Re*(tx*dft/gacobi-ft*dny)-dnx*duxjc-dny*duyjc
        r3=Re*(dpxjc-g1)
        dduxjc=(r1*ts+r2*tm+r3*ty*ty)/gacobi2
        dduyjc=(-r1*ts-r2*tm+r3*tx*tx)/gacobi2
        dduxyjc=(r1*tm-r2*ts-r3*tx*ty)/gacobi2

        r1=-dtx*dvxjc-dty*dvyjc
        r2=-Re*(ty*dft/gacobi+ft*dnx)-dnx*dvxjc-dny*dvyjc
        r3=Re*(dpyjc-g2)
        ddvxjc=(r1*ts+r2*tm+r3*ty*ty)/gacobi2
        ddvyjc=(-r1*ts-r2*tm+r3*tx*tx)/gacobi2
        ddvxyjc=(r1*tm-r2*ts-r3*tx*ty)/gacobi2

        r1=ddfn-dtx*dpxjc-dty*dpyjc
        r2=ddft/gacobi-dgacobi*dft/gacobi2+dgn-
     .     dnx*dpxjc-dny*dpyjc
        r3=0.0d0
        ddpxjc0=(r1*ts+r2*tm+r3*ty*ty)/gacobi2
        ddpyjc0=(-r1*ts-r2*tm+r3*tx*tx)/gacobi2
        r3=1.0d0
        ddpxjc1=(r1*ts+r2*tm+r3*ty*ty)/gacobi2
        ddpyjc1=(-r1*ts-r2*tm+r3*tx*tx)/gacobi2

        ujc(1,m,l)=duxjc
        vjc(1,m,l)=dvxjc
        pjc(1,m,l)=dpxjc
        ujc(3,m,l)=dduxjc
        vjc(3,m,l)=ddvxjc
        pjc(3,m,l)=ddpxjc1
        ujc(5,m,l)=dduxyjc
        vjc(5,m,l)=ddvxyjc
        pjc(5,m,l)=fn
        pjc(7,m,l)=ddpxjc0

        ujc(2,m,l)=duyjc
        vjc(2,m,l)=dvyjc
        pjc(2,m,l)=dpyjc
        ujc(4,m,l)=dduyjc
        vjc(4,m,l)=ddvyjc
        pjc(4,m,l)=ddpyjc1
        ujc(6,m,l)=dduxyjc
        vjc(6,m,l)=ddvxyjc
        pjc(6,m,l)=fn
        pjc(8,m,l)=ddpyjc0
      enddo

      ujc(1,ns,l)=ujc(1,0,l)
      vjc(1,ns,l)=vjc(1,0,l)
      pjc(1,ns,l)=pjc(1,0,l)
      ujc(3,ns,l)=ujc(3,0,l)
      vjc(3,ns,l)=vjc(3,0,l)
      pjc(3,ns,l)=pjc(3,0,l)
      ujc(5,ns,l)=ujc(5,0,l)
      vjc(5,ns,l)=vjc(5,0,l)
      pjc(5,ns,l)=pjc(5,0,l)
      pjc(7,ns,l)=pjc(7,0,l)

      ujc(2,ns,l)=ujc(2,0,l)
      vjc(2,ns,l)=vjc(2,0,l)
      pjc(2,ns,l)=pjc(2,0,l)
      ujc(4,ns,l)=ujc(4,0,l)
      vjc(4,ns,l)=vjc(4,0,l)
      pjc(4,ns,l)=pjc(4,0,l)
      ujc(6,ns,l)=ujc(6,0,l)
      vjc(6,ns,l)=vjc(6,0,l)
      pjc(6,ns,l)=pjc(6,0,l)
      pjc(8,ns,l)=pjc(8,0,l)

      ENDDO

      return
      end


c-----------------------------------------------------------------------
