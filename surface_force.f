c-----------------------------------------------------------------------
c
      subroutine surface_force(fac)
      include 'parameter.inc'
      include 'surface.inc'
      real*8 fac,tx,ty,gacobi,g1,g2,sforce

      sforce=0.0d0

      DO l=1,ms

        do m=0,ns
          tx=taox(m,l)
          ty=taoy(m,l)
          gacobi=dsqrt(tx*tx+ty*ty)
          g1=sforce*(us(m,l)-use(m,l))/(fac*dt)
          g2=sforce*(vs(m,l)-vse(m,l))/(fac*dt)
          gx(m,l)=(g1*tx+g2*ty)/gacobi
          gy(m,l)=(g1*ty-g2*tx)/gacobi
        enddo

      ENDDO

      return
      end


c-----------------------------------------------------------------------
