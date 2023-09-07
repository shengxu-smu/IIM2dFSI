c-----------------------------------------------------------------------
c
      subroutine time_output
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer m1,m2,mm
      real*8 cd(ms),cl(ms)

      m1=0
      m2=ns
      mm=ns/8
      do l=1,ms
        do m=m1,m2,mm
          write(16,100)t,alfa(m),xs(m,l),ys(m,l),xse(m,l),yse(m,l),
     .                 us(m,l),vs(m,l),use(m,l),vse(m,l),
     .                 fx(m,l),fy(m,l)
       enddo
      enddo
100   format(1x,200e16.6)

      do l=1,ms
        cd(l)=0.0d0
        cl(l)=0.0d0
        do m=0,ns-1
          cd(l)=cd(l)-(fx(m,l)*taox(m,l)+fy(m,l)*taoy(m,l))
          cl(l)=cl(l)-(fx(m,l)*taoy(m,l)-fy(m,l)*taox(m,l))
        enddo
        cd(l)=2.0d0*cd(l)*dalfa+2.0d0*area(l)*xsctt(l)
        cl(l)=2.0d0*cl(l)*dalfa+2.0d0*area(l)*ysctt(l)

      enddo
      write(66,100)t,(cd(l),cl(l),l=1,ms),(xsc(l),ysc(l),l=1,ms),
     .             (xsct(l),ysct(l),l=1,ms),(theta(l),thetat(l),l=1,ms)

      return
      end


c-----------------------------------------------------------------------
