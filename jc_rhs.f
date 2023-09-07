c-----------------------------------------------------------------------
c
      subroutine jc_rhs
      include 'parameter.inc'
      include 'surface.inc'
      real*8 foo
      real*8 tmp(2,0:ns),ccc(2,0:ns)

      DO l=1,ms

      do m=0,ns
        tmp(1,m)=fx(m,l)
        tmp(2,m)=fy(m,l)
      enddo
      call cubic_spline(alfa,tmp,ccc,ns,ns,2)
      do m=0,ns-1
        call fdf(alfa(m),alfa(m+1),tmp(1,m),tmp(1,m+1),
     .           ccc(1,m),ccc(1,m+1),alfa(m),foo,dfx(m,l),1)
        call fdf(alfa(m),alfa(m+1),tmp(2,m),tmp(2,m+1),
     .           ccc(2,m),ccc(2,m+1),alfa(m),foo,dfy(m,l),1)
      enddo
      dfx(ns,l)=dfx(0,l)
      dfy(ns,l)=dfy(0,l)
      
      do m=0,ns
        tmp(1,m)=gx(m,l)
        tmp(2,m)=gy(m,l)
      enddo
      call cubic_spline(alfa,tmp,ccc,ns,ns,2)
      do m=0,ns-1
        call fdf(alfa(m),alfa(m+1),tmp(1,m),tmp(1,m+1),
     .           ccc(1,m),ccc(1,m+1),alfa(m),foo,dgx(m,l),1)
        call fdf(alfa(m),alfa(m+1),tmp(2,m),tmp(2,m+1),
     .           ccc(2,m),ccc(2,m+1),alfa(m),foo,dgy(m,l),1)
      enddo
      dgx(ns,l)=dgx(0,l)
      dgy(ns,l)=dgy(0,l)

      do m=0,ns
        tmp(1,m)=dfx(m,l)
        tmp(2,m)=dfy(m,l)
      enddo
      call cubic_spline(alfa,tmp,ccc,ns,ns,2)
      do m=0,ns-1
        call fdf(alfa(m),alfa(m+1),tmp(1,m),tmp(1,m+1),
     .           ccc(1,m),ccc(1,m+1),alfa(m),foo,ddfx(m,l),1)
        call fdf(alfa(m),alfa(m+1),tmp(2,m),tmp(2,m+1),
     .           ccc(2,m),ccc(2,m+1),alfa(m),foo,ddfy(m,l),1)
      enddo
      ddfx(ns,l)=ddfx(0,l)
      ddfy(ns,l)=ddfy(0,l)

      ENDDO

      return
      end


c-----------------------------------------------------------------------
