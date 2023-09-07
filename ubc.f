c-----------------------------------------------------------------------
c
      subroutine ubc(krk,fac)
      include 'parameter.inc'
      include 'field.inc'
      include 'old.inc'
      integer krk
      real*8 fac,conv,visc,uip,uim,vip,vim

c  periodic
      if(ls_ubc.eq.0.and.ln_ubc.eq.0) then
        do i=0,nx
          u(i,0)=u(i,ny-1)
          u(i,ny+1)=u(i,2)
        enddo
      endif

      if(lw_ubc.eq.0.and.le_ubc.eq.0) then
        do j=0,ny+1
          u(0,j)=u(nx-1,j)
          u(nx,j)=u(1,j)
        enddo
      endif

c  dirichlet
      if(lw_ubc.eq.1) then
        do j=0,ny+1
          u(0,j)=-u(1,j)
        enddo
      endif

      if(le_ubc.eq.1) then
        do j=0,ny+1
          u(nx,j)=-u(nx-1,j)
        enddo
      endif

      if(ls_ubc.eq.1) then
        do i=0,nx
          u(i,0)=-u(i,2)
        enddo
      endif

      if(ln_ubc.eq.1) then
        do i=0,nx
          u(i,ny+1)=-u(i,ny-1)
        enddo
      endif

c  newmann
      if(lw_ubc.eq.2) then
        do j=0,ny+1
          u(0,j)=u(1,j)
        enddo
      endif

      if(le_ubc.eq.2) then
        do j=0,ny+1
          u(nx,j)=u(nx-1,j)
        enddo
      endif

      if(ls_ubc.eq.2) then
        do i=0,nx
          u(i,0)=u(i,2)
        enddo
      endif

      if(ln_ubc.eq.2) then
        do i=0,nx
          u(i,ny+1)=u(i,ny-1)
        enddo
      endif

      return
      end


c-----------------------------------------------------------------------
