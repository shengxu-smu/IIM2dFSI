c-----------------------------------------------------------------------
c
      subroutine data_plot
      include 'parameter.inc'
      include 'field.inc'
      integer ih,jh,iref,igap,iend,jgap,jend,ii,jj
      real*8 foo,unorm,vnorm,pnorm
      real*8 au(0:nx,0:ny+1),av(0:nx+1,0:ny),ap(0:nx+1,0:ny+1)

      call streamfunction

      open(unit=29,file='DAT/xe.dat',status='unknown')
      open(unit=39,file='DAT/ye.dat',status='unknown')
      do i=0,nx
        write(29,200)xe(i)
      enddo
      do j=0,ny
        write(39,200)ye(j)
      enddo
      close(29)
      close(39)
      open(unit=49,file='DAT/xc.dat',status='unknown')
      open(unit=59,file='DAT/yc.dat',status='unknown')
      do i=1,nx
        write(49,200)xc(i)
      enddo
      do j=1,ny
        write(59,200)yc(j)
      enddo
      close(49)
      close(59)
200   format(1x,10e16.6)

      ih=(0.0d0-x0)/dx
      jh=(0.0d0-y0)/dy
      open(unit=19,file='DAT/gy.dat',status='unknown')
      open(unit=29,file='DAT/gx.dat',status='unknown')
      do j=1,ny
        write(19,200)yc(j),uxcyc(ih,j),vycxc(ih,j),p(ih,j),o(ih,j)
      enddo
      do i=1,nx
        write(29,200)xc(i),uxcyc(i,jh),vycxc(i,jh),p(i,jh),o(i,jh)
      enddo
      close(19)
      close(29)

      open(unit=38,file='DAT/p.dat',status='unknown')
      open(unit=48,file='DAT/d.dat',status='unknown')
      open(unit=58,file='DAT/uc.dat',status='unknown')
      open(unit=68,file='DAT/vc.dat',status='unknown')
      open(unit=78,file='DAT/wo.dat',status='unknown')
      open(unit=88,file='DAT/ph.dat',status='unknown')
      do j=1,ny
        write(38,400)(p(i,j),i=1,nx)
        write(48,400)(d(i,j),i=1,nx)
        write(58,400)(uxcyc(i,j),i=1,nx)
        write(68,400)(vycxc(i,j),i=1,nx)
        write(78,400)(o(i,j),i=1,nx)
        write(88,400)(-ph(i,j),i=1,nx)
      enddo
      close(38)
      close(48)
      close(58)
      close(68)
      close(78)
      close(88)
400   format(1x,1000e36.18)

      iref=0
      igap=2**2
      iend=igap*(nx-1)+1
      jgap=2**2
      jend=jgap*(ny-1)+1
      if(iref.eq.1.or.iref.eq.-1) then
        open(unit=17,file='DAT/au.dat',status='unknown')
        open(unit=27,file='DAT/av.dat',status='unknown')
        open(unit=37,file='DAT/ap.dat',status='unknown')
      endif
      if(iref.eq.1) then
        do j=1,ny
          write(17,500)(uxcyc(i,j),i=1,nx)
          write(27,500)(vycxc(i,j),i=1,nx)
          write(37,500)(p(i,j),i=1,nx)
        enddo
      endif
      if(iref.eq.-1) then
        do j=1,jend
           do i=1,iend
             if(mod(j-1,jgap).eq.0.and.mod(i-1,igap).eq.0) then
               jj=(j-1)/jgap+1
               ii=(i-1)/igap+1
               read(17,500)au(ii,jj)
               read(27,500)av(ii,jj)
               read(37,500)ap(ii,jj)
             else
               read(17,500)foo
               read(27,500)foo
               read(37,500)foo
             endif
           enddo
         enddo
      endif
      if(iref.eq.1.or.iref.eq.-1) then
        close(17)
        close(27)
        close(37)
      endif

      unorm=0.0d0
      vnorm=0.0d0
      pnorm=0.0d0
      do j=1,ny
        do i=1,nx
          if(iref.eq.-1) then
            uxcyc(i,j)=uxcyc(i,j)-au(i,j)
            vycxc(i,j)=vycxc(i,j)-av(i,j)
            p(i,j)=p(i,j)-ap(i,j)
          endif
          unorm=max(unorm,abs(uxcyc(i,j)))
          vnorm=max(vnorm,abs(vycxc(i,j)))
          pnorm=max(pnorm,abs(p(i,j)))
        enddo
      enddo

      write(*,*)
      write(*,*)'unorm = ',unorm
      write(*,*)'vnorm = ',vnorm
      write(*,*)'pnorm = ',pnorm

      if(iref.eq.-1) then
        open(unit=19,file='DAT/eu.dat',status='unknown')
        open(unit=29,file='DAT/ev.dat',status='unknown')
        open(unit=39,file='DAT/ep.dat',status='unknown')
        do j=1,ny
          write(19,400)(uxcyc(i,j),i=1,nx)
          write(29,400)(vycxc(i,j),i=1,nx)  
          write(39,400)(p(i,j),i=1,nx)
        enddo
        close(19)
        close(29)
        close(39)
      endif
500   format(1x,e36.18)

      return
      end


c-----------------------------------------------------------------------
