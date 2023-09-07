c-----------------------------------------------------------------------
c
      subroutine pbc_fft(pw,pe,ps,pn,krk)
      include 'parameter.inc'
      include 'field.inc'
      integer krk
      real*8 pw(ny),pe(ny),ps(nx),pn(nx)

c  periodic
      do j=1,ny
        pw(j)=0.0d0
        pe(j)=0.0d0
      enddo
      do i=1,nx
        ps(i)=0.0d0
        pn(i)=0.0d0
      enddo

c  dirichlet 
      if(lw_pbc.eq.1) then
        do j=1,ny
          pw(j)=0.0d0
        enddo
      endif

      if(le_pbc.eq.1) then
        do j=1,ny
          pe(j)=0.0d0
        enddo
      endif

      if(ls_pbc.eq.1) then
        do i=1,nx
          ps(i)=0.0d0
        enddo
      endif

      if(ln_pbc.eq.1) then
        do i=1,nx
          pn(i)=0.0d0
        enddo
      endif

c  newmann 
      if(lw_pbc.eq.2) then
        do j=1,ny
          pw(j)=(dx2/Re)*(uxcyc(2,j)-uxcyc(1,j)
     .         +(dx/dy)*(vxeye(0,j)-vxeye(0,j-1)))
c    .         +(dy2/Re)*(uxcyc(1,j+1)-2.0d0*uxcyc(1,j)+uxcyc(1,j-1))
c    .         -dx1*(u(1,j)*u(1,j)-u(0,j)*u(0,j))
c     .         -dy1*(uyexc(1,j)*v(1,j)-uyexc(1,j-1)*v(1,j-1))
        enddo
      endif

      if(le_pbc.eq.2) then
        do j=1,ny
          pe(j)=(dx2/Re)*(uxcyc(nx-1,j)-uxcyc(nx,j)
     .         -(dx/dy)*(vxeye(nx,j)-vxeye(nx,j-1)))
c     .         +(dy2/Re)*(uxcyc(nx,j+1)-2.0d0*uxcyc(nx,j)+uxcyc(nx,j-1))
c     .         -dx1*(u(nx,j)*u(nx,j)-u(nx-1,j)*u(nx-1,j))
c     .         -dy1*(uyexc(nx,j)*v(nx,j)-uyexc(nx,j-1)*v(nx,j-1))
        enddo
      endif

      if(ls_pbc.eq.2) then
        do i=1,nx
          ps(i)=(dy2/Re)*(vycxc(i,2)-vycxc(i,1)
     .         +(dy/dx)*(uyexe(i,0)-uyexe(i-1,0)))
c     .         +(dx2/Re)*(vycxc(i+1,1)-2.0d0*vycxc(i,1)+vycxc(i-1,1))
c     .         -dy1*(v(i,1)*v(i,1)-v(i,0)*v(i,0))
c     .         -dx1*(vycxe(i,1)*u(i,1)-vycxe(i-1,1)*u(i-1,1))
        enddo
      endif

      if(ln_pbc.eq.2) then
        do i=1,nx
          pn(i)=(dy2/Re)*(vycxc(i,ny-1)-vycxc(i,ny)
     .         -(dy/dx)*(uyexe(i,ny)-uyexe(i-1,ny)))
c     .         +(dx2/Re)*(vycxc(i+1,ny)-2.0d0*vycxc(i,ny)+vycxc(i-1,ny))
c     .         -dy1*(v(i,ny)*v(i,ny)-v(i,ny-1)*v(i,ny-1))
c     .         -dx1*(vycxe(i,ny)*u(i,ny)-vycxe(i-1,ny)*u(i-1,ny))
        enddo
      endif

      return
      end


c-----------------------------------------------------------------------
