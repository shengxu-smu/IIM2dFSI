c-----------------------------------------------------------------------
c
      subroutine data_read
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'

      open(unit=11,file='DAT/old.out',form='unformatted',status='old')
      rewind 11
      read(11)nstart,t,t0
      read(11)theta,xsc,ysc,xs,ys,xse,yse
      read(11)x,y,xe,ye,xc,yc
      read(11)us,vs,use,vse,u,v,p,o
      close(11)

      return
      end


c-----------------------------------------------------------------------
c
      subroutine data_write
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'

      open(unit=12,file='DAT/new.out',form='unformatted',
     .     status='unknown')
      rewind 12
      write(12)nend,t,t0
      write(12)theta,xsc,ysc,xs,ys,xse,yse
      write(12)x,y,xe,ye,xc,yc
      write(12)us,vs,use,vse,u,v,p,o
      close(12)

      return
      end


c-----------------------------------------------------------------------
