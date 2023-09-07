c-----------------------------------------------------------------------
c
      subroutine euler_link
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer iu,id,ju,jd
      real*8 alfa0,alfa1,alpha,xsm0,xsm1,ysm0,ysm1,tangent,signx,signy
      real*8 aa(22),bb(22),f(22)
      
      DO l=1,ms

      ncxe(l)=0
      ncxc(l)=0
      ncye(l)=0
      ncyc(l)=0

      do m=0,ns-1

        alfa0=alfa(m)
        alfa1=alfa(m+1)
        xsm0=xs(m,l)
        xsm1=xs(m+1,l)
        ysm0=ys(m,l)
        ysm1=ys(m+1,l)
        aa(1)=(xsm1-xsm0)/dalfa
        bb(1)=(xsm0*alfa1-xsm1*alfa0)/dalfa
        aa(2)=(ysm1-ysm0)/dalfa
        bb(2)=(ysm0*alfa1-ysm1*alfa0)/dalfa
        do n=1,6
          aa(n+2)=(ujc(n,m+1,l)-ujc(n,m,l))/dalfa
          bb(n+2)=(ujc(n,m,l)*alfa1-ujc(n,m+1,l)*alfa0)/dalfa
          aa(n+8)=(vjc(n,m+1,l)-vjc(n,m,l))/dalfa
          bb(n+8)=(vjc(n,m,l)*alfa1-vjc(n,m+1,l)*alfa0)/dalfa
          aa(n+14)=(pjc(n,m+1,l)-pjc(n,m,l))/dalfa
          bb(n+14)=(pjc(n,m,l)*alfa1-pjc(n,m+1,l)*alfa0)/dalfa
        enddo
        aa(21)=(pjc(7,m+1,l)-pjc(7,m,l))/dalfa
        bb(21)=(pjc(7,m,l)*alfa1-pjc(7,m+1,l)*alfa0)/dalfa
        aa(22)=(pjc(8,m+1,l)-pjc(8,m,l))/dalfa
        bb(22)=(pjc(8,m,l)*alfa1-pjc(8,m+1,l)*alfa0)/dalfa

        if(xsm0.gt.xsm1) then
          iu=int((xsm0-x0)/hdx)
          id=int((xsm1-x0)/hdx)+1
          if(iu.ge.id) then
            do i=iu,id,-1
              alpha=(x(i)-bb(1))/aa(1)
              signx=sign(1.0,ysm1-ysm0)
              signy=sign(1.0,xsm0-xsm1)
              if(m.eq.0) then
                tangent=signy+sign(1.0,xs(ns-1,l)-xs(ns,l))
              else
                tangent=signy+sign(1.0,xs(m-1,l)-xsm0)
              endif
              if(x(i).ne.xsm0.or.tangent.ne.0.0) then
                do n=1,22
                  f(n)=aa(n)*alpha+bb(n)
                enddo
                f(2)=ysm0+(ysm1-ysm0)*(x(i)-xsm0)/(xsm1-xsm0)
                j=int((f(2)-y0)/hdy)
                if(f(2).eq.y(j)) j=j-max(0,int(signy))
                if(mod(i,2).eq.0) then
                  ixe(ncxe(l),l)=i/2
                  if(mod(j,2).eq.0) then
                    jexe(ncxe(l),l)=j/2
                    jcxe(ncxe(l),l)=j/2
                  else
                    jexe(ncxe(l),l)=(j-1)/2
                    jcxe(ncxe(l),l)=(j+1)/2
                  endif
                  falfaxe(0,ncxe(l),l)=alpha     
                  falfaxe(1,ncxe(l),l)=f(2)
                  falfaxe(2,ncxe(l),l)=signx
                  falfaxe(3,ncxe(l),l)=signy
                  do n=1,6
                    ujcxe(n,ncxe(l),l)=f(n+2)
                    vjcxe(n,ncxe(l),l)=f(n+8)
                  enddo
                  ncxe(l)=ncxe(l)+1
                else
                  ixc(ncxc(l),l)=(i+1)/2
                  if(mod(j,2).eq.0) then
                    jexc(ncxc(l),l)=j/2
                    jcxc(ncxc(l),l)=j/2
                  else
                    jexc(ncxc(l),l)=(j-1)/2
                    jcxc(ncxc(l),l)=(j+1)/2
                  endif
                  falfaxc(0,ncxc(l),l)=alpha
                  falfaxc(1,ncxc(l),l)=f(2)
                  falfaxc(2,ncxc(l),l)=signx
                  falfaxc(3,ncxc(l),l)=signy
                  do n=1,6
                    ujcxc(n,ncxc(l),l)=f(n+2)
                    vjcxc(n,ncxc(l),l)=f(n+8)
                    pjcxc(n,ncxc(l),l)=f(n+14)
                  enddo
                  pjcxc(7,ncxc(l),l)=f(21)
                  pjcxc(8,ncxc(l),l)=f(22)
                  ncxc(l)=ncxc(l)+1
                endif
              endif
            enddo
          endif
        elseif(xsm1.gt.xsm0) then
          id=int((xsm0-x0)/hdx)+1
          iu=int((xsm1-x0)/hdx)
          if(iu.ge.id) then
            do i=id,iu
              alpha=(x(i)-bb(1))/aa(1)
              signx=sign(1.0,ysm1-ysm0)
              signy=sign(1.0,xsm0-xsm1)
              if(m.eq.ns-1) then
                tangent=signy+sign(1.0,xs(0,l)-xs(1,l))
              else
                tangent=signy+sign(1.0,xsm1-xs(m+2,l))
              endif
              if(x(i).ne.xsm1.or.tangent.ne.0.0) then
                do n=1,22
                  f(n)=aa(n)*alpha+bb(n)
                enddo
                f(2)=ysm0+(ysm1-ysm0)*(x(i)-xsm0)/(xsm1-xsm0)
                j=int((f(2)-y0)/hdy)
                if(f(2).eq.y(j)) j=j-max(0,int(signy))
                if(mod(i,2).eq.0) then
                  ixe(ncxe(l),l)=i/2
                  if(mod(j,2).eq.0) then
                    jexe(ncxe(l),l)=j/2
                    jcxe(ncxe(l),l)=j/2
                  else
                    jexe(ncxe(l),l)=(j-1)/2
                    jcxe(ncxe(l),l)=(j+1)/2
                  endif
                  falfaxe(0,ncxe(l),l)=alpha
                  falfaxe(1,ncxe(l),l)=f(2)
                  falfaxe(2,ncxe(l),l)=signx
                  falfaxe(3,ncxe(l),l)=signy
                  do n=1,6
                    ujcxe(n,ncxe(l),l)=f(n+2)
                    vjcxe(n,ncxe(l),l)=f(n+8)
                  enddo
                  ncxe(l)=ncxe(l)+1
                else
                  ixc(ncxc(l),l)=(i+1)/2
                  if(mod(j,2).eq.0) then
                    jexc(ncxc(l),l)=j/2
                    jcxc(ncxc(l),l)=j/2
                  else
                    jexc(ncxc(l),l)=(j-1)/2
                    jcxc(ncxc(l),l)=(j+1)/2
                  endif
                  falfaxc(0,ncxc(l),l)=alpha
                  falfaxc(1,ncxc(l),l)=f(2)
                  falfaxc(2,ncxc(l),l)=signx
                  falfaxc(3,ncxc(l),l)=signy
                  do n=1,6
                    ujcxc(n,ncxc(l),l)=f(n+2)
                    vjcxc(n,ncxc(l),l)=f(n+8)
                    pjcxc(n,ncxc(l),l)=f(n+14)
                  enddo
                  pjcxc(7,ncxc(l),l)=f(21)
                  pjcxc(8,ncxc(l),l)=f(22)
                  ncxc(l)=ncxc(l)+1
                endif
              endif
            enddo
          endif
        endif

        if(ysm0.gt.ysm1) then
          ju=int((ysm0-y0)/hdy)
          jd=int((ysm1-y0)/hdy)+1
          if(ju.ge.jd) then
            do j=ju,jd,-1
              alpha=(y(j)-bb(2))/aa(2)
              signx=sign(1.0,ysm1-ysm0)
              signy=sign(1.0,xsm0-xsm1)
              if(m.eq.0) then
                tangent=signx+sign(1.0,ys(ns,l)-ys(ns-1,l))
              else
                tangent=signx+sign(1.0,ysm0-ys(m-1,l))
              endif
              if(y(j).ne.ysm0.or.tangent.ne.0.0) then
                do n=1,22
                  f(n)=aa(n)*alpha+bb(n)
                enddo
                f(1)=xsm0+(xsm1-xsm0)*(y(j)-ysm0)/(ysm1-ysm0)
                i=int((f(1)-x0)/hdx)
                if(f(1).eq.x(i)) i=i-max(0,int(signx))
                if(mod(j,2).eq.0) then
                  jye(ncye(l),l)=j/2
                  if(mod(i,2).eq.0) then
                    ieye(ncye(l),l)=i/2
                    icye(ncye(l),l)=i/2
                  else
                    ieye(ncye(l),l)=(i-1)/2
                    icye(ncye(l),l)=(i+1)/2
                  endif
                  falfaye(0,ncye(l),l)=alpha
                  falfaye(1,ncye(l),l)=f(1)
                  falfaye(2,ncye(l),l)=signx
                  falfaye(3,ncye(l),l)=signy
                  do n=1,6
                    ujcye(n,ncye(l),l)=f(n+2)
                    vjcye(n,ncye(l),l)=f(n+8)
                  enddo
                  ncye(l)=ncye(l)+1
                else
                  jyc(ncyc(l),l)=(j+1)/2
                  if(mod(i,2).eq.0) then
                    ieyc(ncyc(l),l)=i/2
                    icyc(ncyc(l),l)=i/2
                  else
                    ieyc(ncyc(l),l)=(i-1)/2
                    icyc(ncyc(l),l)=(i+1)/2
                  endif
                  falfayc(0,ncyc(l),l)=alpha
                  falfayc(1,ncyc(l),l)=f(1)
                  falfayc(2,ncyc(l),l)=signx
                  falfayc(3,ncyc(l),l)=signy
                  do n=1,6
                    ujcyc(n,ncyc(l),l)=f(n+2)
                    vjcyc(n,ncyc(l),l)=f(n+8)
                    pjcyc(n,ncyc(l),l)=f(n+14)
                  enddo
                  pjcyc(7,ncyc(l),l)=f(21)
                  pjcyc(8,ncyc(l),l)=f(22)
                  ncyc(l)=ncyc(l)+1
                endif
              endif
            enddo
          endif
        elseif(ysm1.gt.ysm0) then
          jd=int((ysm0-y0)/hdy)+1
          ju=int((ysm1-y0)/hdy)
          if(ju.ge.jd) then
            do j=jd,ju
              alpha=(y(j)-bb(2))/aa(2)
              signx=sign(1.0,ysm1-ysm0)
              signy=sign(1.0,xsm0-xsm1)
              if(m.eq.ns-1) then
                tangent=signx+sign(1.0,ys(1,l)-ys(0,l))
              else
                tangent=signx+sign(1.0,ys(m+2,l)-ysm1)
              endif
              if(y(j).ne.ysm1.or.tangent.ne.0.0) then
                do n=1,22
                  f(n)=aa(n)*alpha+bb(n)
                enddo
                f(1)=xsm0+(xsm1-xsm0)*(y(j)-ysm0)/(ysm1-ysm0)
                i=int((f(1)-x0)/hdx)
                if(f(1).eq.x(i)) i=i-max(0,int(signx))
                if(mod(j,2).eq.0) then
                  jye(ncye(l),l)=j/2
                  if(mod(i,2).eq.0) then
                    ieye(ncye(l),l)=i/2
                    icye(ncye(l),l)=i/2
                  else
                    ieye(ncye(l),l)=(i-1)/2
                    icye(ncye(l),l)=(i+1)/2
                  endif
                  falfaye(0,ncye(l),l)=alpha
                  falfaye(1,ncye(l),l)=f(1)
                  falfaye(2,ncye(l),l)=signx
                  falfaye(3,ncye(l),l)=signy
                  do n=1,6
                    ujcye(n,ncye(l),l)=f(n+2)
                    vjcye(n,ncye(l),l)=f(n+8)
                  enddo
                  ncye(l)=ncye(l)+1
                else
                  jyc(ncyc(l),l)=(j+1)/2
                  if(mod(i,2).eq.0) then
                    ieyc(ncyc(l),l)=i/2
                    icyc(ncyc(l),l)=i/2
                  else
                    ieyc(ncyc(l),l)=(i-1)/2
                    icyc(ncyc(l),l)=(i+1)/2
                  endif
                  falfayc(0,ncyc(l),l)=alpha
                  falfayc(1,ncyc(l),l)=f(1)
                  falfayc(2,ncyc(l),l)=signx
                  falfayc(3,ncyc(l),l)=signy
                  do n=1,6
                    ujcyc(n,ncyc(l),l)=f(n+2)
                    vjcyc(n,ncyc(l),l)=f(n+8)
                    pjcyc(n,ncyc(l),l)=f(n+14)
                  enddo
                  pjcyc(7,ncyc(l),l)=f(21)
                  pjcyc(8,ncyc(l),l)=f(22)
                  ncyc(l)=ncyc(l)+1
                endif
              endif
            enddo
          endif
        endif

      enddo

      ncxe(l)=ncxe(l)-1
      ncxc(l)=ncxc(l)-1
      ncye(l)=ncye(l)-1
      ncyc(l)=ncyc(l)-1

      ENDDO

      return
      end


c-----------------------------------------------------------------------
