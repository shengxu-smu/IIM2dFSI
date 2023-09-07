c-----------------------------------------------------------------------
c
      subroutine poisson_iterate(lw_bc,le_bc,ls_bc,ln_bc,
     .                           nx1,nx2,ny1,ny2,dx,dy,phi,zeta,krk)
      integer krk
      real*8 pi,error
      parameter(pi=
     .3.1415926535897932384626433832795028841971693993751058209749446d0)
      parameter (error=1.0d-12) 
      integer lw_bc,le_bc,ls_bc,ln_bc,nx1,nx2,ny1,ny2,Ix,Jy,i,j,iter
      real*8 dx,dy,beta,xi,omega,dold,dnew
      real*8 phi(nx1:nx2,ny1:ny2),zeta(nx1:nx2,ny1:ny2)

      Ix=nx2-nx1-1
      Jy=ny2-ny1-1
      beta=dx/dy
      xi=(dcos(pi/float(Ix))+
     .    beta*beta*dcos(pi/float(Jy)))/(1.0d0+beta*beta)
      omega=2.0d0/(1.0d0+dsqrt(1.0d0-xi*xi))

      if(lw_bc.eq.0.and.le_bc.eq.0) then
        do j=ny1,ny2
          phi(nx1,j)=phi(nx2-1,j)
          phi(nx2,j)=phi(nx1+1,j)
        enddo
      endif
      if(ln_bc.eq.0.and.ls_bc.eq.0) then
        do i=nx1,nx2
          phi(i,ny1)=phi(i,ny2-1)
          phi(i,ny2)=phi(i,ny1+1)
        enddo
      endif

      dold=0.0
      do i=nx1+1,nx2-1
        do j=ny1+1,ny2-1
          dold=dold+abs(phi(i,j))
        enddo
      enddo

      iter=0
 5    iter=iter+1
      dnew=0.0

      phi(nx1+1,ny1+1)=phi(nx1+1,ny1+1)+
     .                 (omega/(2.0d0*(1.0d0+beta*beta)))*
     .                 (phi(nx1+2,ny1+1)+phi(nx1,ny1+1)+
     .                 beta*beta*(phi(nx1+1,ny1+2)+phi(nx1+1,ny1))-
     .                 2.0d0*(1.0d0+beta*beta)*phi(nx1+1,ny1+1)-
     .                 dx*dx*zeta(nx1+1,ny1+1))
      if(lw_bc.eq.0.and.le_bc.eq.0.and.ln_bc.eq.0.and.ls_bc.eq.0) then
        phi(nx1+1,ny1+1)=0.0d0
      endif
      if(lw_bc.eq.0.and.le_bc.eq.0) then
        phi(nx2,ny1+1)=phi(nx1+1,ny1+1)
      endif
      if(ln_bc.eq.0.and.ls_bc.eq.0) then
        phi(nx1+1,ny2)=phi(nx1+1,ny1+1)
      endif
      dnew=dnew+abs(phi(nx1+1,ny1+1))

      do i=nx1+2,nx2-1
        phi(i,ny1+1)=phi(i,ny1+1)+(omega/(2.0d0*(1.0d0+beta*beta)))*
     .               (phi(i+1,ny1+1)+phi(i-1,ny1+1)+
     .               beta*beta*(phi(i,ny1+2)+phi(i,ny1))-
     .               2.0d0*(1.0d0+beta*beta)*phi(i,ny1+1)-
     .               dx*dx*zeta(i,ny1+1))
        if(ln_bc.eq.0.and.ls_bc.eq.0) then
          phi(i,ny2)=phi(i,ny1+1)
        endif
        dnew=dnew+abs(phi(i,ny1+1))
      enddo

      do j=ny1+2,ny2-1
        phi(nx1+1,j)=phi(nx1+1,j)+(omega/(2.0d0*(1.0d0+beta*beta)))*
     .               (phi(nx1+2,j)+phi(nx1,j)+
     .               beta*beta*(phi(nx1+1,j+1)+phi(nx1+1,j-1))-
     .               2.0d0*(1.0d0+beta*beta)*phi(nx1+1,j)-
     .               dx*dx*zeta(nx1+1,j))
        if(lw_bc.eq.0.and.le_bc.eq.0) then
          phi(nx2,j)=phi(nx1+1,j)
        endif
        dnew=dnew+abs(phi(nx1+1,j))
      enddo

      do j=ny1+2,ny2-1
        do i=nx1+2,nx2-1
          phi(i,j)=phi(i,j)+(omega/(2.0d0*(1.0d0+beta*beta)))*
     .             (phi(i+1,j)+phi(i-1,j)+
     .             beta*beta*(phi(i,j+1)+phi(i,j-1))-
     .             2.0d0*(1.0d0+beta*beta)*phi(i,j)-
     .             dx*dx*zeta(i,j))
          dnew=dnew+abs(phi(i,j))
        enddo
      enddo

      if(lw_bc.eq.0.and.le_bc.eq.0) then
        do j=ny1+1,ny2-1
          phi(nx1,j)=phi(nx2-1,j)
        enddo
      endif
      if(ln_bc.eq.0.and.ls_bc.eq.0) then
        do i=nx1+1,nx2-1
          phi(i,ny1)=phi(i,ny2-1)
        enddo
      endif

      if(lw_bc.eq.2.or.le_bc.eq.2.or.ln_bc.eq.2.or.ls_bc.eq.2) then
        call pbc_iterate(nx1,nx2,ny1,ny2,phi,krk)
      endif

      if(lw_bc.eq.3.or.le_bc.eq.3.or.ln_bc.eq.3.or.ls_bc.eq.3) then
        call pbc_iterate(nx1,nx2,ny1,ny2,phi,krk)
      endif

      if(abs(dnew-dold).ge.error) then
        dold=dnew
        goto 5
      endif
      write(*,*)'   !!! iter = ',iter

      return
      end


c-----------------------------------------------------------------------
