c-----------------------------------------------------------------------
c
      subroutine singular_force
      include 'parameter.inc'
      include 'surface.inc'
      real*8 foo,alpha,tx,ty,gacobi2,gacobie2,s,se,curv,tension,f1,f2
      real*8 coord(2,0:ns),ccc(2,0:ns),taoxe(0:ns),taoye(0:ns)

      DO l=1,ms

        do m=0,ns
          fx(m,l)=0.0d0
          fy(m,l)=0.0d0
        enddo

        do m=0,ns
          coord(1,m)=xse(m,l)
          coord(2,m)=yse(m,l)
        enddo
        call cubic_spline(alfa,coord,ccc,ns,ns,2)
        do m=0,ns-1
          call fdf(alfa(m),alfa(m+1),coord(1,m),coord(1,m+1),
     .             ccc(1,m),ccc(1,m+1),alfa(m),foo,taoxe(m),1)
          call fdf(alfa(m),alfa(m+1),coord(2,m),coord(2,m+1),
     .             ccc(2,m),ccc(2,m+1),alfa(m),foo,taoye(m),1)
        enddo
        taoxe(ns)=taoxe(0)
        taoye(ns)=taoye(0)

        do m=0,ns
          tx=taox(m,l)
          ty=taoy(m,l)
          gacobi2=tx*tx+ty*ty
          gacobie2=taoxe(m)*taoxe(m)+taoye(m)*taoye(m)
          tension=sk0(l)*(dsqrt(gacobi2/gacobie2)-1.0d0)
          coord(1,m)=tension*tx
          coord(2,m)=tension*ty
        enddo
        call cubic_spline(alfa,coord,ccc,ns,ns,2)
        do m=0,ns-1
          call fdf(alfa(m),alfa(m+1),coord(1,m),coord(1,m+1),
     .             ccc(1,m),ccc(1,m+1),alfa(m),foo,fx(m,l),1)
          call fdf(alfa(m),alfa(m+1),coord(2,m),coord(2,m+1),
     .             ccc(2,m),ccc(2,m+1),alfa(m),foo,fy(m,l),1)
        enddo
        fx(ns,l)=fx(0,l)
        fy(ns,l)=fy(0,l)

        do m=0,ns
          tx=taox(m,l)
          ty=taoy(m,l)
          gacobi2=tx*tx+ty*ty
          if(ifxy.eq.1) then
            s=(xs(m,l)-xsc(l))**2.0d0+(ys(m,l)-ysc(l))**2.0d0
            se=(xse(m,l)-xsc(l))**2.0d0+(yse(m,l)-ysc(l))**2.0d0
            f1=fx(m,l)-sk1(l)*(us(m,l)-use(m,l))
     .                -sk2(l)*(xs(m,l)-xse(m,l))
     .                -sk3(l)*xs(m,l)*(1.0d0-dsqrt(se/s))
            f2=fy(m,l)-sk1(l)*(vs(m,l)-vse(m,l))
     .                -sk2(l)*(ys(m,l)-yse(m,l))
     .                -sk3(l)*ys(m,l)*(1.0d0-dsqrt(se/s))
            fx(m,l)=(f1*tx+f2*ty)/gacobi2
            fy(m,l)=(f1*ty-f2*tx)/gacobi2
          else
            curv=-(dtaoy(m,l)*tx-dtaox(m,l)*ty)/gacobi2**1.5d0
            fx(m,l)=0.0d0
            fy(m,l)=0.2d0*curv
          endif
        enddo

      ENDDO

      return
      end


c-----------------------------------------------------------------------
