c-----------------------------------------------------------------------
c
      subroutine body_force(fac)
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer nn
      real*8 fac,bforce,xx,yy,ang,r,rc

      bforce=0.0d0

      do j=0,ny+1
        do i=0,nx
          gu(i,j)=0.0d0
        enddo
      enddo
      do j=0,ny
        do i=0,nx+1
          gv(i,j)=0.0d0
        enddo
      enddo

      DO l=1,ms

      do j=1,ny
        do i=1,nx-1
          xx=xe(i)-xsc(l)
          yy=yc(j)-ysc(l)
          ang=datan2(yy,xx)+pi-theta(l)
          if(ang.lt.0.0d0) ang=ang+2.0d0*pi
          r=dsqrt(xx*xx+yy*yy)
          nn=int(ang/dar)
          rc=ra(nn,l)
          if(r.le.rc) then
            gu(i,j)=bforce*(xsct(l)-thetat(l)*yy-u(i,j))/(fac*dt)
          endif
        enddo
      enddo

      do j=1,ny-1
        do i=1,nx
          xx=xc(i)-xsc(l)
          yy=ye(j)-ysc(l)
          ang=datan2(yy,xx)+pi-theta(l)
          if(ang.lt.0.0d0) ang=ang+2.0d0*pi
          r=dsqrt(xx*xx+yy*yy)
          nn=int(ang/dar)
          rc=ra(nn,l)
          if(r.le.rc) then
            gv(i,j)=bforce*(ysct(l)+thetat(l)*xx-v(i,j))/(fac*dt)
          endif
        enddo
      enddo

      ENDDO

      return
      end


c-----------------------------------------------------------------------
