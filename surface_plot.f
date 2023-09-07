c-----------------------------------------------------------------------
c
      subroutine surface_plot
      include 'parameter.inc'
      include 'surface.inc'
      character*1 fnext
      character*16 fname

      DO l=1,ms

      write(unit=fnext,fmt='(i1)') l
      fname='DAT/surface'//fnext//'.dat'
      open(unit=69,file=fname,status='unknown')
      do m=0,ns
        write(69,300)alfa(m),xs(m,l),ys(m,l),xse(m,l),yse(m,l),
     .               taox(m,l),taoy(m,l),dtaox(m,l),dtaoy(m,l),
     .               us(m,l),vs(m,l),use(m,l),vse(m,l),
     .               fx(m,l),fy(m,l),dfx(m,l),dfy(m,l),
     .               ddfx(m,l),ddfy(m,l),
     .               (ujc(n,m,l),vjc(n,m,l),pjc(n,m,l),n=1,6)
      enddo
      close(69)

      ENDDO

300   format(1x,40e16.6)

      return
      end


c-----------------------------------------------------------------------
