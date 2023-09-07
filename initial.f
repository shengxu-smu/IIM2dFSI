c-----------------------------------------------------------------------
c
      subroutine initial
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer ni,nj,ip
      real*8 tc

      t=0.0d0
      t0=0.0d0
      ip=1

      if(ip.eq.1) area(1)=a(1)*b(1)+pi*b(1)*b(1)/4.0d0
      if(ip.eq.3) area(1)=pi*a(1)*b(1)

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
 
      xsctt(1)=1.25d0*(-0.64d0*dcos(0.8d0*t0))*dcos(pi/3.0d0)
      ysctt(1)=1.25d0*(-0.64d0*dcos(0.8d0*t0))*dsin(pi/3.0d0)

      ENDDO

      call vrffti(ns,wsave)
      call surface_parametrization

      DO l=1,ms
        do m=0,ns
          us(m,l)=0.0d0
          vs(m,l)=0.0d0
          use(m,l)=xsct(l)-thetat(l)*(yse(m,l)-ysc(l))
          vse(m,l)=ysct(l)+thetat(l)*(xse(m,l)-xsc(l))
        enddo
      ENDDO

      do j=0,ny+1
        do i=0,nx
          u(i,j)=0.0d0
          gu(i,j)=0.0d0
        enddo
      enddo

      do j=0,ny
        do i=0,nx+1
          v(i,j)=0.0d0
          gv(i,j)=0.0d0
        enddo
      enddo

      do j=0,ny+1
        do i=0,nx+1
          p(i,j)=0.0d0
        enddo
      enddo

      do j=0,ny+1
        do i=0,nx+1
          d(i,j)=0.0d0
        enddo
      enddo

      if(lw_pbc.eq.0.and.le_pbc.eq.0.and.
     .   ls_pbc.eq.0.and.ln_pbc.eq.0) then
        ni=nx-1
        nj=ny-1
        call vrffti(ni,wsavei)
        call vrffti(nj,wsavej)
      endif

      if(lw_pbc.eq.1.and.le_pbc.eq.1.and.
     .   ls_pbc.eq.1.and.ln_pbc.eq.1) then
        ni=nx-2
        nj=ny-2
        call vsinti(ni,wsavei)
        call vsinti(nj,wsavej)
      endif

      if(lw_pbc.eq.2.and.le_pbc.eq.2.and.
     .   ls_pbc.eq.2.and.ln_pbc.eq.2) then
        ni=nx
        nj=ny
        call vcosti(ni,wsavei)
        call vcosti(nj,wsavej)
      endif

      if(lw_pbc.eq.1.and.le_pbc.eq.2.and.
     .   ls_pbc.eq.2.and.ln_pbc.eq.2) then
        ni=nx-1
        nj=ny
        call vsinqi(ni,wsavei)
        call vcosti(nj,wsavej)
      endif      

      if(lw_pbc.eq.1.and.le_pbc.eq.2.and.
     .   ls_pbc.eq.0.and.ln_pbc.eq.0) then
        ni=nx-1
        nj=ny-1
        call vsinqi(ni,wsavei)
        call vrffti(nj,wsavej)
      endif

      return
      end


c-----------------------------------------------------------------------
