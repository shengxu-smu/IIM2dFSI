c field.inc
      real*8 x(-2:2*nx+2),y(-2:2*ny+2)
      real*8 xe(-1:nx+1),ye(-1:ny+1),xc(0:nx+1),yc(0:ny+1)
      real*8 wsavei(3*nx+15),wsavej(3*ny+15)

      real*8 u(0:nx,0:ny+1),v(0:nx+1,0:ny),p(0:nx+1,0:ny+1)
      real*8 gu(0:nx,0:ny+1),gv(0:nx+1,0:ny)
      real*8 d(0:nx+1,0:ny+1),o(0:nx+1,0:ny+1),ph(0:nx+1,0:ny+1)
      real*8 uyexe(0:nx,0:ny),uyexc(nx,0:ny),uxcyc(nx,0:ny+1)
      real*8 vxeye(0:nx,0:ny),vycxe(0:nx,ny),vycxc(0:nx+1,ny)

      common/grid/x,y,xe,ye,xc,yc
      common/field/u,v,p,d,o,ph,gu,gv,
     .             uyexe,uyexc,uxcyc,vxeye,vycxe,vycxc
      common/fftfield/wsavei,wsavej
