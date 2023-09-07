c-----------------------------------------------------------------------
c
      subroutine field_interpolate
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'

      do j=0,ny
        do i=0,nx
          uyexe(i,j)=0.5d0*(u(i,j)+u(i,j+1))
          vxeye(i,j)=0.5d0*(v(i,j)+v(i+1,j))
        enddo
      enddo
      if(isingular.eq.1) then
        do l=1,ms
          do j=0,ncye(l)
            vxeye(ivxeye(j,l),jye(j,l))=vxeye(ivxeye(j,l),jye(j,l))+
     .                                  vcixeye(j,l)
          enddo
        enddo
      endif

      do j=0,ny+1
        do i=1,nx
          uxcyc(i,j)=0.5d0*(u(i-1,j)+u(i,j))
        enddo
      enddo
      if(isingular.eq.1) then
        do l=1,ms
          do j=0,ncyc(l)
            uxcyc(iuxcyc(j,l),jyc(j,l))=uxcyc(iuxcyc(j,l),jyc(j,l))+
     .                                  ucixcyc(j,l)
          enddo
        enddo
      endif

      do j=0,ny
        do i=1,nx
          uyexc(i,j)=0.5d0*(uxcyc(i,j)+uxcyc(i,j+1))
        enddo
      enddo

      do j=1,ny
        do i=0,nx+1
          vycxc(i,j)=0.5d0*(v(i,j-1)+v(i,j))
        enddo
      enddo

      do j=1,ny
        do i=0,nx
          vycxe(i,j)=0.5d0*(vxeye(i,j-1)+vxeye(i,j))
        enddo
      enddo

      if(isingular.eq.1) then
        do l=1,ms
          do i=0,ncxe(l)
            uyexe(ixe(i,l),juyexe(i,l))=uyexe(ixe(i,l),juyexe(i,l))+
     .                                  uciyexe(i,l)
            vycxe(ixe(i,l),jvycxe(i,l))=vycxe(ixe(i,l),jvycxe(i,l))+
     .                                  vciycxe(i,l)
          enddo
          do i=0,ncxc(l)
            uyexc(ixc(i,l),juyexc(i,l))=uyexc(ixc(i,l),juyexc(i,l))+
     .                                  uciyexc(i,l)
            vycxc(ixc(i,l),jvycxc(i,l))=vycxc(ixc(i,l),jvycxc(i,l))+
     .                                  vciycxc(i,l)
          enddo
        enddo
      endif

      return
      end


c-----------------------------------------------------------------------
