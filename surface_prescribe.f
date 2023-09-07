c-----------------------------------------------------------------------  
c
      subroutine surface_prescribe(krk,fac)
      include 'parameter.inc'
      include 'surface.inc'
      integer krk
      real*8 fac,tc

      if(krk.eq.1) t0=t0+0.5d0*dt
      if(krk.eq.3) t0=t0+0.5d0*dt

      tc=1.0d0

      DO l=1,ms

      theta(l)=pi-0.25d0*pi*(1.0d0-dsin(0.8d0*t0)*
     .                      (1.0d0-dexp(-t0/tc)))
      xsc(l)=1.25d0*(dcos(0.8d0*t0)+1.0d0)*dcos(pi/3.0d0)
      ysc(l)=1.25d0*(dcos(0.8d0*t0)+1.0d0)*dsin(pi/3.0d0)

      thetat(l)=0.25d0*pi*(0.8d0*dcos(0.8d0*t0))*(1.0d0-dexp(-t0/tc))+
     .          0.25d0*pi*(dsin(0.8d0*t0))*(dexp(-t0/tc)/tc)
      xsct(l)=1.25d0*(-0.8d0*dsin(0.8d0*t0))*dcos(pi/3.0d0)
      ysct(l)=1.25d0*(-0.8d0*dsin(0.8d0*t0))*dsin(pi/3.0d0)

      xsctt(l)=1.25d0*(-0.64d0*dcos(0.8d0*t0))*dcos(pi/3.0d0)
      ysctt(l)=1.25d0*(-0.64d0*dcos(0.8d0*t0))*dsin(pi/3.0d0)

      do m=0,ns-1
        xse(m,l)=xsc(l)+xs0(m,l)*dcos(theta(l))-ys0(m,l)*dsin(theta(l))
        yse(m,l)=ysc(l)+xs0(m,l)*dsin(theta(l))+ys0(m,l)*dcos(theta(l))
        use(m,l)=xsct(l)-thetat(l)*(yse(m,l)-ysc(l))
        vse(m,l)=ysct(l)+thetat(l)*(xse(m,l)-xsc(l))
      enddo
      xse(ns,l)=xse(0,l)
      yse(ns,l)=yse(0,l)
      use(ns,l)=use(0,l)
      vse(ns,l)=vse(0,l)

      ENDDO

      return
      end


c-----------------------------------------------------------------------
