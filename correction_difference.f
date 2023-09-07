c-----------------------------------------------------------------------
c
      subroutine correction_difference
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer ie,ie1,ic,ic1,je,je1,jc,jc1
      real*8 sx,sy,hp,hm,signx,signy
      real*8 ui,vi,djc,ddjc


      DO l=1,ms

      do i=0,ncxe(l)
        signy=falfaxe(3,i,l)
        sy=falfaxe(1,i,l)
        ie=ixe(i,l)
        je=jexe(i,l)
        je1=je+1
        jc=jcxe(i,l)
        jc1=jc+1
        if(jc.eq.je) then
          hp=yc(jc1)-sy
          hm=sy-ye(je)
          ui=(hp*uyexe(ie,je)+hm*u(ie,jc1))/hdy-
     .       ujcxe(2,i,l)*hp*hm/hdy-signy*dy1*hp*hm*ujcxe(4,i,l)*
     .       ((1.0d0+signy)*hp+(1.0d0-signy)*hm)
          vi=(hp*vxeye(ie,je)+hm*vycxe(ie,jc1))/hdy-
     .       vjcxe(2,i,l)*hp*hm/hdy-signy*dy1*hp*hm*vjcxe(4,i,l)*
     .       ((1.0d0+signy)*hp+(1.0d0-signy)*hm)

          djc=ui*vjcxe(2,i,l)+vi*ujcxe(2,i,l)
          ddjc=ui*vjcxe(4,i,l)+vi*ujcxe(4,i,l)+2.0d0*
     .         (falfaxe(4,i,l)*falfaxe(6,i,l)-
     .          falfaxe(5,i,l)*falfaxe(7,i,l))

          ucdy(i,l)=-dy1*(djc*(ye(je)-sy)+
     .              0.5d0*ddjc*(ye(je)-sy)**2.0d0)
          ucdyy(1,i,l)=-dy2*(ujcxe(2,i,l)*(yc(jc1)-sy)+
     .                 0.5d0*ujcxe(4,i,l)*(yc(jc1)-sy)**2.0d0)
          ucdyy(2,i,l)=dy2*(ujcxe(2,i,l)*(yc(jc)-sy)+
     .                0.5d0*ujcxe(4,i,l)*(yc(jc)-sy)**2.0d0)
        else
          hp=ye(je1)-sy
          hm=sy-yc(jc)
          ui=(hp*u(ie,jc)+hm*uyexe(ie,je1))/hdy-
     .       ujcxe(2,i,l)*hp*hm/hdy-signy*dy1*hp*hm*ujcxe(4,i,l)*
     .       ((1.0d0+signy)*hp+(1.0d0-signy)*hm)
          vi=(hp*vycxe(ie,jc)+hm*vxeye(ie,je1))/hdy-
     .       vjcxe(2,i,l)*hp*hm/hdy-signy*dy1*hp*hm*vjcxe(4,i,l)*
     .       ((1.0d0+signy)*hp+(1.0d0-signy)*hm)

          djc=ui*vjcxe(2,i,l)+vi*ujcxe(2,i,l)
          ddjc=ui*vjcxe(4,i,l)+vi*ujcxe(4,i,l)+2.0d0*
     .         (falfaxe(4,i,l)*falfaxe(6,i,l)-
     .          falfaxe(5,i,l)*falfaxe(7,i,l))

          ucdy(i,l)=-dy1*(djc*(ye(je1)-sy)+
     .             0.5d0*ddjc*(ye(je1)-sy)**2.0d0)
          ucdyy(1,i,l)=-dy2*(ujcxe(2,i,l)*(yc(jc1)-sy)+
     .                 0.5d0*ujcxe(4,i,l)*(yc(jc1)-sy)**2.0d0)
          ucdyy(2,i,l)=dy2*(ujcxe(2,i,l)*(yc(jc)-sy)+
     .                0.5d0*ujcxe(4,i,l)*(yc(jc)-sy)**2.0d0)
        endif
        falfaxe(13,i,l)=ui
        falfaxe(14,i,l)=vi
      enddo

      do i=0,ncxc(l)
        signy=falfaxc(3,i,l)
        sy=falfaxc(1,i,l)
        ic=ixc(i,l)
        je=jexc(i,l)
        je1=je+1
        jc=jcxc(i,l)
        jc1=jc+1
        if(jc.eq.je) then
          hp=yc(jc1)-sy
          hm=sy-ye(je)
          ui=(hp*uyexc(ic,je)+hm*uxcyc(ic,jc1))/hdy-
     .       ujcxc(2,i,l)*hp*hm/hdy-signy*dy1*hp*hm*ujcxc(4,i,l)*
     .       ((1.0d0+signy)*hp+(1.0d0-signy)*hm)
          vi=(hp*v(ic,je)+hm*vycxc(ic,jc1))/hdy-
     .       vjcxc(2,i,l)*hp*hm/hdy-signy*dy1*hp*hm*vjcxc(4,i,l)*
     .       ((1.0d0+signy)*hp+(1.0d0-signy)*hm)

          djc=2.0d0*vi*vjcxc(2,i,l)
          ddjc=2.0d0*vi*vjcxc(4,i,l)+2.0d0*
     .         (falfaxc(10,i,l)*falfaxc(10,i,l)-
     .          falfaxc(11,i,l)*falfaxc(11,i,l))

          vcvdy(i,l)=-dy1*(djc*(yc(jc1)-sy)+
     .              0.5d0*ddjc*(yc(jc1)-sy)**2.0d0)
          vcpdy(i,l)=-dy1*(pjcxc(6,i,l)+
     .                     pjcxc(2,i,l)*(yc(jc1)-sy)+
     .               0.5d0*pjcxc(4,i,l)*(yc(jc1)-sy)**2.0d0)
          vcdyy(1,i,l)=-dy2*(vjcxc(2,i,l)*(ye(je1)-sy)+
     .                 0.5d0*vjcxc(4,i,l)*(ye(je1)-sy)**2.0d0)
          vcdyy(2,i,l)=dy2*(vjcxc(2,i,l)*(ye(je)-sy)+
     .                0.5d0*vjcxc(4,i,l)*(ye(je)-sy)**2.0d0)
          pcdyy(1,i,l)=-dy2*(pjcxc(6,i,l)+
     .                       pjcxc(2,i,l)*(yc(jc1)-sy)+
     .                 0.5d0*pjcxc(4,i,l)*(yc(jc1)-sy)**2.0d0)
          pcdyy(2,i,l)=dy2*(pjcxc(6,i,l)+
     .                      pjcxc(2,i,l)*(yc(jc)-sy)+
     .                0.5d0*pjcxc(4,i,l)*(yc(jc)-sy)**2.0d0)
        else
          hp=ye(je1)-sy
          hm=sy-yc(jc)
          ui=(hp*uxcyc(ic,jc)+hm*uyexc(ic,je1))/hdy-
     .       ujcxc(2,i,l)*hp*hm/hdy-signy*dy1*hp*hm*ujcxc(4,i,l)*
     .       ((1.0d0+signy)*hp+(1.0d0-signy)*hm)
          vi=(hp*vycxc(ic,jc)+hm*v(ic,je1))/hdy-
     .       vjcxc(2,i,l)*hp*hm/hdy-signy*dy1*hp*hm*vjcxc(4,i,l)*
     .       ((1.0d0+signy)*hp+(1.0d0-signy)*hm)

          djc=2.0d0*vi*vjcxc(2,i,l)
          ddjc=2.0d0*vi*vjcxc(4,i,l)+2.0d0*
     .         (falfaxc(10,i,l)*falfaxc(10,i,l)-
     .          falfaxc(11,i,l)*falfaxc(11,i,l))

          vcvdy(i,l)=-dy1*(djc*(yc(jc)-sy)+
     .              0.5d0*ddjc*(yc(jc)-sy)**2.0d0)
          vcpdy(i,l)=-dy1*(pjcxc(6,i,l)+
     .                     pjcxc(2,i,l)*(yc(jc)-sy)+
     .               0.5d0*pjcxc(4,i,l)*(yc(jc)-sy)**2.0d0)
          vcdyy(1,i,l)=-dy2*(vjcxc(2,i,l)*(ye(je1)-sy)+
     .                 0.5d0*vjcxc(4,i,l)*(ye(je1)-sy)**2.0d0)
          vcdyy(2,i,l)=dy2*(vjcxc(2,i,l)*(ye(je)-sy)+
     .                0.5d0*vjcxc(4,i,l)*(ye(je)-sy)**2.0d0)
          pcdyy(1,i,l)=-dy2*(pjcxc(6,i,l)+
     .                       pjcxc(2,i,l)*(yc(jc1)-sy)+
     .                 0.5d0*pjcxc(4,i,l)*(yc(jc1)-sy)**2.0d0)
          pcdyy(2,i,l)=dy2*(pjcxc(6,i,l)+
     .                      pjcxc(2,i,l)*(yc(jc)-sy)+
     .                0.5d0*pjcxc(4,i,l)*(yc(jc)-sy)**2.0d0)
        endif
        falfaxc(13,i,l)=ui
        falfaxc(14,i,l)=vi
      enddo

      do j=0,ncye(l)
        signx=falfaye(2,j,l)
        sx=falfaye(1,j,l)
        je=jye(j,l)
        ie=ieye(j,l)
        ie1=ie+1
        ic=icye(j,l)
        ic1=ic+1
        if(ic.eq.ie) then
          hp=xc(ic1)-sx
          hm=sx-xe(ie)
          ui=(hp*uyexe(ie,je)+hm*uyexc(ic1,je))/hdx-
     .       ujcye(1,j,l)*hp*hm/hdx-signx*dx1*hp*hm*ujcye(3,j,l)*
     .       ((1.0d0+signx)*hp+(1.0d0-signx)*hm)
          vi=(hp*vxeye(ie,je)+hm*v(ic1,je))/hdx-
     .       vjcye(1,j,l)*hp*hm/hdx-signx*dx1*hp*hm*vjcye(3,j,l)*
     .       ((1.0d0+signx)*hp+(1.0d0-signx)*hm)

          djc=ui*vjcye(1,j,l)+vi*ujcye(1,j,l)
          ddjc=ui*vjcye(3,j,l)+vi*ujcye(3,j,l)+2.0d0*
     .         (falfaye(4,j,l)*falfaye(6,j,l)-
     .          falfaye(5,j,l)*falfaye(7,j,l))

          vcdx(j,l)=-dx1*(djc*(xe(ie)-sx)+
     .             0.5d0*ddjc*(xe(ie)-sx)**2.0d0)
          vcdxx(1,j,l)=-dx2*(vjcye(1,j,l)*(xc(ic1)-sx)+
     .                 0.5d0*vjcye(3,j,l)*(xc(ic1)-sx)**2.0d0)
          vcdxx(2,j,l)=dx2*(vjcye(1,j,l)*(xc(ic)-sx)+
     .                0.5d0*vjcye(3,j,l)*(xc(ic)-sx)**2.0d0)
        else
          hp=xe(ie1)-sx
          hm=sx-xc(ic)
          ui=(hp*uyexc(ic,je)+hm*uyexe(ie1,je))/hdx-
     .       ujcye(1,j,l)*hp*hm/hdx-signx*dx1*hp*hm*ujcye(3,j,l)*
     .       ((1.0d0+signx)*hp+(1.0d0-signx)*hm)
          vi=(hp*v(ic,je)+hm*vxeye(ie1,je))/hdx-
     .       vjcye(1,j,l)*hp*hm/hdx-signx*dx1*hp*hm*vjcye(3,j,l)*
     .       ((1.0d0+signx)*hp+(1.0d0-signx)*hm)

          djc=ui*vjcye(1,j,l)+vi*ujcye(1,j,l)
          ddjc=ui*vjcye(3,j,l)+vi*ujcye(3,j,l)+2.0d0*
     .         (falfaye(4,j,l)*falfaye(6,j,l)-
     .          falfaye(5,j,l)*falfaye(7,j,l))

          vcdx(j,l)=-dx1*(djc*(xe(ie1)-sx)+
     .             0.5d0*ddjc*(xe(ie1)-sx)**2.0d0)
          vcdxx(1,j,l)=-dx2*(vjcye(1,j,l)*(xc(ic1)-sx)+
     .                 0.5d0*vjcye(3,j,l)*(xc(ic1)-sx)**2.0d0)
          vcdxx(2,j,l)=dx2*(vjcye(1,j,l)*(xc(ic)-sx)+
     .                0.5d0*vjcye(3,j,l)*(xc(ic)-sx)**2.0d0)
        endif
        falfaye(13,j,l)=ui
        falfaye(14,j,l)=vi
      enddo

      do j=0,ncyc(l)
        signx=falfayc(2,j,l)
        sx=falfayc(1,j,l)
        jc=jyc(j,l)
        ie=ieyc(j,l)
        ie1=ie+1
        ic=icyc(j,l)
        ic1=ic+1
        if(ic.eq.ie) then
          hp=xc(ic1)-sx
          hm=sx-xe(ie)
          ui=(hp*u(ie,jc)+hm*uxcyc(ic1,jc))/hdx-
     .       ujcyc(1,j,l)*hp*hm/hdx-signx*dx1*hp*hm*ujcyc(3,j,l)*
     .       ((1.0d0+signx)*hp+(1.0d0-signx)*hm)
          vi=(hp*vycxe(ie,jc)+hm*vycxc(ic1,jc))/hdx-
     .       vjcyc(1,j,l)*hp*hm/hdx-signx*dx1*hp*hm*vjcyc(3,j,l)*
     .       ((1.0d0+signx)*hp+(1.0d0-signx)*hm)

          djc=2.0d0*ui*ujcyc(1,j,l)
          ddjc=2.0d0*ui*ujcyc(3,j,l)+2.0d0*
     .         (falfayc(4,j,l)*falfayc(4,j,l)-
     .          falfayc(5,j,l)*falfayc(5,j,l))

          ucudx(j,l)=-dx1*(djc*(xc(ic1)-sx)+
     .              0.5d0*ddjc*(xc(ic1)-sx)**2.0d0)
          ucpdx(j,l)=-dx1*(pjcyc(5,j,l)+
     .                     pjcyc(1,j,l)*(xc(ic1)-sx)+
     .               0.5d0*pjcyc(3,j,l)*(xc(ic1)-sx)**2.0d0)
          ucdxx(1,j,l)=-dx2*(ujcyc(1,j,l)*(xe(ie1)-sx)+
     .                 0.5d0*ujcyc(3,j,l)*(xe(ie1)-sx)**2.0d0)
          ucdxx(2,j,l)=dx2*(ujcyc(1,j,l)*(xe(ie)-sx)+
     .                0.5d0*ujcyc(3,j,l)*(xe(ie)-sx)**2.0d0)
          pcdxx(1,j,l)=-dx2*(pjcyc(5,j,l)+
     .                       pjcyc(1,j,l)*(xc(ic1)-sx)+
     .                 0.5d0*pjcyc(3,j,l)*(xc(ic1)-sx)**2.0d0)
          pcdxx(2,j,l)=dx2*(pjcyc(5,j,l)+
     .                      pjcyc(1,j,l)*(xc(ic)-sx)+
     .                0.5d0*pjcyc(3,j,l)*(xc(ic)-sx)**2.0d0)
        else
          hp=xe(ie1)-sx
          hm=sx-xc(ic)
          ui=(hp*uxcyc(ic,jc)+hm*u(ie1,jc))/hdx-
     .       ujcyc(1,j,l)*hp*hm/hdx-signx*dx1*hp*hm*ujcyc(3,j,l)*
     .       ((1.0d0+signx)*hp+(1.0d0-signx)*hm)
          vi=(hp*vycxc(ic,jc)+hm*vycxe(ie1,jc))/hdx-
     .       vjcyc(1,j,l)*hp*hm/hdx-signx*dx1*hp*hm*vjcyc(3,j,l)*
     .       ((1.0d0+signx)*hp+(1.0d0-signx)*hm)

          djc=2.0d0*ui*ujcyc(1,j,l)
          ddjc=2.0d0*ui*ujcyc(3,j,l)+2.0d0*
     .         (falfayc(4,j,l)*falfayc(4,j,l)-
     .          falfayc(5,j,l)*falfayc(5,j,l))

          ucudx(j,l)=-dx1*(djc*(xc(ic)-sx)+
     .              0.5d0*ddjc*(xc(ic)-sx)**2.0d0)
          ucpdx(j,l)=-dx1*(pjcyc(5,j,l)+
     .                     pjcyc(1,j,l)*(xc(ic)-sx)+
     .               0.5d0*pjcyc(3,j,l)*(xc(ic)-sx)**2.0d0)
          ucdxx(1,j,l)=-dx2*(ujcyc(1,j,l)*(xe(ie1)-sx)+
     .                 0.5d0*ujcyc(3,j,l)*(xe(ie1)-sx)**2.0d0)
          ucdxx(2,j,l)=dx2*(ujcyc(1,j,l)*(xe(ie)-sx)+
     .                0.5d0*ujcyc(3,j,l)*(xe(ie)-sx)**2.0d0)
          pcdxx(1,j,l)=-dx2*(pjcyc(5,j,l)+
     .                       pjcyc(1,j,l)*(xc(ic1)-sx)+
     .                 0.5d0*pjcyc(3,j,l)*(xc(ic1)-sx)**2.0d0)
          pcdxx(2,j,l)=dx2*(pjcyc(5,j,l)+
     .                      pjcyc(1,j,l)*(xc(ic)-sx)+
     .                0.5d0*pjcyc(3,j,l)*(xc(ic)-sx)**2.0d0)
        endif
        falfayc(13,j,l)=ui
        falfayc(14,j,l)=vi
      enddo

      ENDDO

      return
      end


c-----------------------------------------------------------------------
