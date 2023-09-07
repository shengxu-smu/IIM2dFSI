c-----------------------------------------------------------------------
c
      subroutine lagrangian_interpolate
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer ncx,ncy,nc,nn,ip,nsamex,nsamey,nsamexy,ifilter
      real*8 bl,bu,gacobi
      real*8 ax(0:2*nx),bx(0:2*nx),abx(0:4*nx)
      real*8 fax(2,0:2*nx),fbx(2,0:2*nx),fabx(2,0:4*nx)
      real*8 ay(0:2*ny),by(0:2*ny),aby(0:4*ny)
      real*8 fay(2,0:2*ny),fby(2,0:2*ny),faby(2,0:4*ny)
      real*8 abxy(0:4*(nx+ny)),fabxy(2,0:4*(nx+ny)),ccxy(2,0:4*(nx+ny))
      real*8 ffm(2),ffmp(2),ccm(2),ccmp(2),foo(2),dfoo(2)
      real*8 r(2,ns),w(2,ns)

      DO l=1,ms

      do i=0,ncxe(l)
         ax(i)=falfaxe(0,i,l)
         fax(1,i)=falfaxe(13,i,l)
         fax(2,i)=falfaxe(14,i,l)
      enddo
      do i=0,ncxc(l)
         bx(i)=falfaxc(0,i,l)
         fbx(1,i)=falfaxc(13,i,l)
         fbx(2,i)=falfaxc(14,i,l)
      enddo
      call order_combine(ax,ncxe(l),2*nx,bx,ncxc(l),2*nx,abx,
     .                   fax,fbx,fabx,2,nsamex)

      do j=0,ncye(l)
         ay(j)=falfaye(0,j,l)
         fay(1,j)=falfaye(13,j,l)
         fay(2,j)=falfaye(14,j,l)
      enddo
      do j=0,ncyc(l)
         by(j)=falfayc(0,j,l)
         fby(1,j)=falfayc(13,j,l)
         fby(2,j)=falfayc(14,j,l)
      enddo
      call order_combine(ay,ncye(l),2*ny,by,ncyc(l),2*ny,aby,
     .                   fay,fby,faby,2,nsamey)

      ncx=ncxe(l)+ncxc(l)+1-nsamex
      ncy=ncye(l)+ncyc(l)+1-nsamey
      call order_combine(abx,ncx,4*nx,aby,ncy,4*ny,abxy,
     .                   fabx,faby,fabxy,2,nsamexy)

      nn=4*(nx+ny)
      nc=ncx+ncy+1-nsamexy      
      if(abs(abxy(nc)-abxy(0)-2.0d0*pi).gt.1.0d-5) then
        nc=nc+1
        abxy(nc)=abxy(0)+2.0d0*pi
        fabxy(1,nc)=fabxy(1,0)
        fabxy(2,nc)=fabxy(2,0)
      endif

      call cubic_spline(abxy,fabxy,ccxy,nc,nn,2)

      if(abxy(0).gt.0.0d0) then
        ip=-1
      else
        ip=0
      endif
      do m=0,ns-1
5       continue
        if(ip.eq.-1) then
          bl=abxy(nc-1)-2.0d0*pi
          if(abs(bl).lt.1.0d-5) bl=0.0d0
          bu=abxy(0)
          do n=1,2
            ffm(n)=fabxy(n,nc-1)
            ffmp(n)=fabxy(n,nc)
            ccm(n)=ccxy(n,nc-1)
            ccmp(n)=ccxy(n,nc)
          enddo
        else
          bl=abxy(ip)
          bu=abxy(ip+1)
          do n=1,2
            ffm(n)=fabxy(n,ip)
            ffmp(n)=fabxy(n,ip+1)
            ccm(n)=ccxy(n,ip)
            ccmp(n)=ccxy(n,ip+1)
          enddo          
        endif
        if(alfa(m).ge.bl.and.alfa(m).lt.bu) then
          call fdf(bl,bu,ffm,ffmp,
     .             ccm,ccmp,alfa(m),foo,dfoo,2)
          us(m,l)=foo(1)
          vs(m,l)=foo(2)
        else
          ip=ip+1
          goto 5              
        endif
      enddo
      us(ns,l)=us(0,l)
      vs(ns,l)=vs(0,l)

      ifilter=1
      if(ifilter.eq.1) then
        do m=0,ns-1
          r(1,m+1)=taox(m,l)*us(m,l)+taoy(m,l)*vs(m,l)
          r(2,m+1)=taoy(m,l)*us(m,l)-taox(m,l)*vs(m,l)
        enddo

        call vrfftf(2,ns,r,w,2,wsave)
        r(2,1)=0.0d0
        do m=64,ns
          r(1,m)=0.0d0
          r(2,m)=0.0d0
        enddo
        call vrfftb(2,ns,r,w,2,wsave)

        do m=0,ns-1
          gacobi=taox(m,l)*taox(m,l)+taoy(m,l)*taoy(m,l)
          us(m,l)=(taox(m,l)*r(1,m+1)+taoy(m,l)*r(2,m+1))/gacobi
          vs(m,l)=(taoy(m,l)*r(1,m+1)-taox(m,l)*r(2,m+1))/gacobi
        enddo
        us(ns,l)=us(0,l)
        vs(ns,l)=vs(0,l)
      endif

      ENDDO

      return
      end


c-----------------------------------------------------------------------
