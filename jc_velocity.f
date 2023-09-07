c-----------------------------------------------------------------------
c
      subroutine jc_velocity
      include 'parameter.inc'
      include 'surface.inc'
      real*8 signx,signy

      DO l=1,ms

      do i=0,ncxe(l)
        signx=falfaxe(2,i,l)
        signy=falfaxe(3,i,l)

        ujcxe(1,i,l)=signx*ujcxe(1,i,l)
        vjcxe(1,i,l)=signx*vjcxe(1,i,l)
        ujcxe(3,i,l)=signx*ujcxe(3,i,l)
        vjcxe(3,i,l)=signx*vjcxe(3,i,l)
        ujcxe(5,i,l)=signx*ujcxe(5,i,l)
        vjcxe(5,i,l)=signx*vjcxe(5,i,l)

        ujcxe(2,i,l)=signy*ujcxe(2,i,l)
        vjcxe(2,i,l)=signy*vjcxe(2,i,l)
        ujcxe(4,i,l)=signy*ujcxe(4,i,l)
        vjcxe(4,i,l)=signy*vjcxe(4,i,l)
        ujcxe(6,i,l)=signy*ujcxe(6,i,l)
        vjcxe(6,i,l)=signy*vjcxe(6,i,l)
      enddo

      do i=0,ncxc(l)
        signx=falfaxc(2,i,l)
        signy=falfaxc(3,i,l)

        ujcxc(1,i,l)=signx*ujcxc(1,i,l)
        vjcxc(1,i,l)=signx*vjcxc(1,i,l)
        ujcxc(3,i,l)=signx*ujcxc(3,i,l)
        vjcxc(3,i,l)=signx*vjcxc(3,i,l)
        ujcxc(5,i,l)=signx*ujcxc(5,i,l)
        vjcxc(5,i,l)=signx*vjcxc(5,i,l)

        ujcxc(2,i,l)=signy*ujcxc(2,i,l)
        vjcxc(2,i,l)=signy*vjcxc(2,i,l)
        ujcxc(4,i,l)=signy*ujcxc(4,i,l)
        vjcxc(4,i,l)=signy*vjcxc(4,i,l)
        ujcxc(6,i,l)=signy*ujcxc(6,i,l)
        vjcxc(6,i,l)=signy*vjcxc(6,i,l)
      enddo

      do j=0,ncye(l)
        signx=falfaye(2,j,l)
        signy=falfaye(3,j,l)

        ujcye(1,j,l)=signx*ujcye(1,j,l)
        vjcye(1,j,l)=signx*vjcye(1,j,l)
        ujcye(3,j,l)=signx*ujcye(3,j,l)
        vjcye(3,j,l)=signx*vjcye(3,j,l)
        ujcye(5,j,l)=signx*ujcye(5,j,l)
        vjcye(5,j,l)=signx*vjcye(5,j,l)

        ujcye(2,j,l)=signy*ujcye(2,j,l)
        vjcye(2,j,l)=signy*vjcye(2,j,l)
        ujcye(4,j,l)=signy*ujcye(4,j,l)
        vjcye(4,j,l)=signy*vjcye(4,j,l)
        ujcye(6,j,l)=signy*ujcye(6,j,l)
        vjcye(6,j,l)=signy*vjcye(6,j,l)
      enddo

      do j=0,ncyc(l)
        signx=falfayc(2,j,l)
        signy=falfayc(3,j,l)

        ujcyc(1,j,l)=signx*ujcyc(1,j,l)
        vjcyc(1,j,l)=signx*vjcyc(1,j,l)
        ujcyc(3,j,l)=signx*ujcyc(3,j,l)
        vjcyc(3,j,l)=signx*vjcyc(3,j,l)
        ujcyc(5,j,l)=signx*ujcyc(5,j,l)
        vjcyc(5,j,l)=signx*vjcyc(5,j,l)

        ujcyc(2,j,l)=signy*ujcyc(2,j,l)
        vjcyc(2,j,l)=signy*vjcyc(2,j,l)
        ujcyc(4,j,l)=signy*ujcyc(4,j,l)
        vjcyc(4,j,l)=signy*vjcyc(4,j,l)
        ujcyc(6,j,l)=signy*ujcyc(5,j,l)
        vjcyc(6,j,l)=signy*vjcyc(6,j,l)
      enddo

      ENDDO

      return
      end


c-----------------------------------------------------------------------
