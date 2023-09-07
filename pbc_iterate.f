c-----------------------------------------------------------------------
c
      subroutine pbc_iterate(nx1,nx2,ny1,ny2,phi,krk)
      include 'parameter.inc'
      include 'field.inc'
      integer nx1,nx2,ny1,ny2,krk
      real*8 phi(nx1:nx2,ny1:ny2)
      real*8 tmp

      if(krk.eq.1) tmp=t-dt
      if(krk.eq.2) tmp=t-0.5d0*dt
      if(krk.eq.3) tmp=t-0.5d0*dt
      if(krk.eq.4) tmp=t

c  periodic
      if(lw_pbc.eq.0.and.le_pbc.eq.0) then
        do j=ny1,ny2
          phi(nx1,j)=phi(nx2-1,j)
          phi(nx2,j)=phi(nx1+1,j)
        enddo
      endif

      if(ln_pbc.eq.0.and.ls_pbc.eq.0) then
        do i=nx1,nx2
          phi(i,ny1)=phi(i,ny2-1)
          phi(i,ny2)=phi(i,ny1+1)
        enddo
      endif

c  dirichlet
      if(lw_pbc.eq.1) then
        do j=ny1,ny2
          phi(nx1,j)=0.0d0
        enddo
      endif

      if(le_pbc.eq.1) then
        do j=ny1,ny2
          phi(nx2,j)=0.0d0
        enddo
      endif

      if(ls_pbc.eq.1) then
        do i=nx1,nx2
          phi(i,ny1)=0.0d0
        enddo
      endif

      if(ln_pbc.eq.1) then
        do i=nx1,nx2
          phi(i,ny2)=0.0d0
        enddo
      endif

c  newmann
      if(lw_pbc.eq.2) then
        do j=ny1,ny2
          phi(nx1,j)=phi(nx1+1,j)
        enddo
      endif

      if(le_pbc.eq.2) then
        do j=ny1,ny2          
          phi(nx2,j)=phi(nx2-1,j)+(1.0d0/(Re*dx))*
     .               (u(nx2-1,j)-(dx/dy)*(v(nx2-1,j)-v(nx2-1,j-1))-
     .               2.0d0*u(nx2-1,j)+u(nx2-2,j))
        enddo
      endif

      if(ln_pbc.eq.2) then
        do i=nx1,nx2
          phi(i,ny2)=phi(i,ny2-1)
        enddo
      endif

      if(ls_pbc.eq.2) then
        do i=nx1,nx2
          phi(i,ny1)=phi(i,ny1+1)
        enddo
      endif
        
      return
      end


c-----------------------------------------------------------------------
