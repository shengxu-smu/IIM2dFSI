c-----------------------------------------------------------------------  
c
      subroutine surface_move(krk,fac)
      include 'parameter.inc'
      include 'surface.inc'
      include 'old.inc'
      integer krk,ifilter,motion
      real*8 fac
      real*8 r(2,ns),w(2,ns)

      call lagrangian_interpolate

      DO l=1,ms

      do m=0,ns-1
        usrk(krk,m,l)=us(m,l)
        vsrk(krk,m,l)=vs(m,l)
        if(krk.eq.4) then
          xs(m,l)=xsn(m,l)+(dt/6.0d0)*(usrk(1,m,l)+
     .            2.0d0*(usrk(2,m,l)+usrk(3,m,l))+usrk(4,m,l))
          ys(m,l)=ysn(m,l)+(dt/6.0d0)*(vsrk(1,m,l)+
     .            2.0d0*(vsrk(2,m,l)+vsrk(3,m,l))+vsrk(4,m,l))
        else
          xs(m,l)=xsn(m,l)+fac*dt*us(m,l)
          ys(m,l)=ysn(m,l)+fac*dt*vs(m,l)
        endif
        r(1,m+1)=xs(m,l)
        r(2,m+1)=ys(m,l)
      enddo
      xs(ns,l)=xs(0,l)
      ys(ns,l)=ys(0,l)

      ifilter=1
      if(ifilter.eq.1) then
        call vrfftf(2,ns,r,w,2,wsave)      
        do m=32,ns
          r(1,m)=0.0d0
          r(2,m)=0.0d0
        enddo
        call vrfftb(2,ns,r,w,2,wsave)

        do m=0,ns-1
          xs(m,l)=r(1,m+1)
          ys(m,l)=r(2,m+1)
        enddo
        xs(ns,l)=xs(0,l)
        ys(ns,l)=ys(0,l)
      endif

      ENDDO

      motion=1
      if(motion.eq.1) call surface_prescribe(krk,fac)

      return
      end


c-----------------------------------------------------------------------
