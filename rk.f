c-----------------------------------------------------------------------
c
      subroutine rk
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      include 'old.inc'
      integer iflag,krk
      real*8 fac

c  krk=1:
      if(isingular.eq.1) then
        call correction_interpolate
        call field_interpolate
        call correction_strain
        call uv_strain
        call dudv_surface
        call jc_pressure
        call correction_difference
        call body_force(0.5d0)
      else
        call field_interpolate
        call uv_strain
        call body_force(0.5d0)
      endif
      
      call old_save

      krk=1
      fac=0.5d0

      if(move.eq.1) then
        call surface_move(krk,fac)
      endif

      call pressure(krk,fac)

      call rhs(krk)
 1    call pgrad(krk)

      do j=1,ny
        do i=1,nx-1
          u(i,j)=un(i,j)+fac*dt*urk(i,j,1)
        enddo
      enddo
      do j=1,ny-1
        do i=1,nx
          v(i,j)=vn(i,j)+fac*dt*vrk(i,j,1)
        enddo
      enddo 

      call ubc(krk,fac)
      call vbc(krk,fac)

      if(isingular.eq.1) then
        call surface_property
        call singular_force
        call surface_force(0.5d0)
        call jc_rhs
        call jc_firstsecond
        call euler_link
        call jc_velocity
        if(itemporal.eq.1) then
          call time_correct(krk,fac)
        endif
      endif

      if(ichorin.eq.1) then
        call chorin_iterate(krk,fac,iflag)
        if(iflag.eq.1) goto 1
      endif

      return

c  krk=2:
      if(isingular.eq.1) then
        call correction_interpolate
        call field_interpolate
        call correction_strain
        call uv_strain
        call dudv_surface
        call jc_pressure
        call correction_difference
        call body_force(0.5d0)
      else
        call field_interpolate
        call uv_strain
        call body_force(0.5d0)
      endif

      krk=2
      fac=0.5d0

      if(move.eq.1) then
        call surface_move(krk,fac)
      endif

      call pressure(krk,fac)

      call rhs(krk)
 2    call pgrad(krk)

      do j=1,ny
        do i=1,nx-1
          u(i,j)=un(i,j)+fac*dt*urk(i,j,2)
        enddo
      enddo
      do j=1,ny-1
        do i=1,nx
          v(i,j)=vn(i,j)+fac*dt*vrk(i,j,2)
        enddo
      enddo

      call ubc(krk,fac)
      call vbc(krk,fac)

      if(isingular.eq.1) then
        call surface_property
        call singular_force
        call surface_force(1.0d0)
        call jc_rhs
        call jc_firstsecond
        call euler_link
        call jc_velocity
        if(itemporal.eq.1) then
          call time_correct(krk,fac)
        endif
      endif

      if(ichorin.eq.1) then
        call chorin_iterate(krk,fac,iflag)
        if(iflag.eq.1) goto 2
      endif

c  krk=3:
      if(isingular.eq.1) then
        call correction_interpolate
        call field_interpolate
        call correction_strain
        call uv_strain
        call dudv_surface
        call jc_pressure
        call correction_difference
        call body_force(1.0d0)
      else
        call field_interpolate
        call uv_strain
        call body_force(1.0d0)
      endif

      krk=3 
      fac=1.0d0

      if(move.eq.1) then
        call surface_move(krk,fac)
      endif

      call pressure(krk,fac)

      call rhs(krk)
 3    call pgrad(krk)

      do j=1,ny
        do i=1,nx-1
          u(i,j)=un(i,j)+fac*dt*urk(i,j,3)
        enddo
      enddo
      do j=1,ny-1
        do i=1,nx
          v(i,j)=vn(i,j)+fac*dt*vrk(i,j,3)
        enddo
      enddo

      call ubc(krk,fac)
      call vbc(krk,fac)

      if(isingular.eq.1) then
        call surface_property
        call singular_force
        call surface_force(1.0d0)
        call jc_rhs
        call jc_firstsecond
        call euler_link
        call jc_velocity
        if(itemporal.eq.1) then
          call time_correct(krk,fac)
        endif
      endif

      if(ichorin.eq.1) then
        call chorin_iterate(krk,fac,iflag)
        if(iflag.eq.1) goto 3
      endif

c  krk=4:
      if(isingular.eq.1) then
        call correction_interpolate
        call field_interpolate
        call correction_strain
        call uv_strain
        call dudv_surface
        call jc_pressure
        call correction_difference
        call body_force(0.5d0)
      else
        call field_interpolate
        call uv_strain
        call body_force(0.5d0)
      endif

      krk=4
      fac=1.0d0

      if(move.eq.1) then
        call surface_move(krk,fac)
      endif

      call pressure(krk,fac)

      call rhs(krk)
 4    call pgrad(krk)

      do j=1,ny
        do i=1,nx-1
          u(i,j)=un(i,j)+(dt/6.0d0)*(urk(i,j,1)+
     .                   2.0d0*(urk(i,j,2)+urk(i,j,3))+urk(i,j,4))
        enddo
      enddo
      do j=1,ny-1
        do i=1,nx
          v(i,j)=vn(i,j)+(dt/6.0d0)*(vrk(i,j,1)+
     .                   2.0d0*(vrk(i,j,2)+vrk(i,j,3))+vrk(i,j,4))
        enddo
      enddo

      call ubc(krk,fac)
      call vbc(krk,fac)

      if(isingular.eq.1) then
        call surface_property
        call singular_force
        call surface_force(0.5d0)
        call jc_rhs
        call jc_firstsecond
        call euler_link
        call jc_velocity
        if(itemporal.eq.1) then
          call time_correct(krk,fac)
        endif
      endif

      if(ichorin.eq.1) then
        call chorin_iterate(krk,fac,iflag)
        if(iflag.eq.1) goto 4
      endif

      return
      end


c-----------------------------------------------------------------------
