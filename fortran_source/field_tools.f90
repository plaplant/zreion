module field_tools
  ! Imports for C linking
  use iso_c_binding

  ! Cosmology
  use cosmo_tools

  ! Domain tools
  use domain_tools

  ! Particle tools
  use particle_tools

  ! Default
  implicit none


contains


  subroutine calc_density_velocity_fields(z, delta1, density, velocity, &
       dep_scheme)
    ! Subroutine arguments
    real(c_double), intent(in)    :: z
    real(c_double), intent(in)    :: delta1(:,:,:)
    real(c_double), intent(inout) :: density(:,:,:),velocity(:,:,:,:)
    integer(c_int), intent(in)    :: dep_scheme
    ! Local variables
    integer(4)            :: i,j,k,n
    real(8)               :: x_dm_gr,mass,ngridr,ncell
    real(8)               :: d,davg,dstd,dmax,dmin,vel_unit
    real(8), dimension(3) :: v,vavg,vstd,vmax,vmin

    ! Unpack parameters
    x_dm_gr = real(Ngrid, kind=8) / real(Ndm, kind=8)
    mass    = x_dm_gr**3
    ngridr  = real(Ngrid, kind=8)
    ncell   = ngridr**3

    ! Check that grid sizes match expectations
    if (size(density, dim=1) /= Ngrid .or. &
         size(density, dim=2) /= Ngrid .or. &
         size(density, dim=3) /= Ngrid) then
       write(*,*) "Bad size for density array"
       return
    endif

    if (size(velocity, dim=1) /= 3 .or. &
         size(velocity, dim=2) /= Ngrid .or. &
         size(velocity, dim=3) /= Ngrid .or. &
         size(velocity, dim=4) /= Ngrid) then
       write(*,*) "Bad size for velocity array"
       return
    endif

    ! See if we need to (re)initialize domain
    if (.not. allocated(domain) .or. &
         size(domain, dim=2) /= Ndomain) then
       call init_domain
    endif

    ! See if we need to (re)initialize dm
    if (.not. allocated(dm) .or. &
         size(dm, dim=2) /= ndm) then
       call init_dm
    endif

    ! Calculate cosmology
    call calc_cosmo_vars(z)
    vel_unit   = x_unit / time_unit

    ! Compute 2LPT
    call calc_2lpt(z, delta1)

    ! Displace particles
    call calc_dm

    ! Construct linked list
    call link_dm

    ! Initialize fields
    !$omp parallel do     &
    !$omp default(shared) &
    !$omp private(k)
    do k=1,Ngrid
       density(:, :, k)     = 0
       velocity(:, :, :, k) = 0
    enddo
    !$omp end parallel do

    select case(dep_scheme)
    case(0)  ! NGP
       !$omp parallel        &
       !$omp default(shared) &
       !$omp private(i,n)
       !$omp do schedule(dynamic, 1)
       do n=1,Ndomain
          call set_domain(i, Ndomain, domain)
          call ngp_mass_assignment(indx_dom(:, :, i))
          call end_domain(i, Ndomain, domain)
       enddo
       !$omp end do
       !$omp do schedule(dynamic, 1)
       do n=1,Ndomain
          call set_domain(i, Ndomain, domain)
          call ngp_vel_assignment(indx_dom(:, :, i))
          call end_domain(i, Ndomain, domain)
       enddo
       !$omp end do
       !$omp end parallel
    case(1)  ! CIC
       !$omp parallel        &
       !$omp default(shared) &
       !$omp private(i,n)
       !$omp do schedule(dynamic, 1)
       do n=1,Ndomain
          call set_domain(i, Ndomain, domain)
          call cic_mass_assignment(indx_dom(:, :, i))
          call end_domain(i, Ndomain, domain)
       enddo
       !$omp end do
       !$omp do schedule(dynamic, 1)
       do n=1,Ndomain
          call set_domain(i, Ndomain, domain)
          call cic_vel_assignment(indx_dom(:, :, i))
          call end_domain(i, Ndomain, domain)
       enddo
       !$omp end do
       !$omp end parallel
    case(2)  ! TSC
       !$omp parallel        &
       !$omp default(shared) &
       !$omp private(i,n)
       !$omp do schedule(dynamic, 1)
       do n=1,Ndomain
          call set_domain(i, Ndomain, domain)
          call tsc_mass_assignment(indx_dom(:, :, i))
          call end_domain(i, Ndomain, domain)
       enddo
       !$omp end do
       !$omp do schedule(dynamic, 1)
       do n=1,Ndomain
          call set_domain(i, Ndomain, domain)
          call tsc_vel_assignment(indx_dom(:, :, i))
          call end_domain(i, Ndomain, domain)
       enddo
       !$omp end do
       !$omp end parallel
    end select


    ! Calculate array statistics
    ! Initialize
    davg = 0
    dstd = 0
    dmax = -huge(dmax)
    dmin =  huge(dmin)

    vavg = 0
    vstd = 0
    vmax = -huge(vmax)
    vmin =  huge(vmin)

    ! Compute in parallel
    !$omp parallel do                      &
    !$omp default(shared)                  &
    !$omp private(i,j,k,d,v)               &
    !$omp reduction(+:davg,dstd,vavg,vstd) &
    !$omp reduction(max:dmax,vmax)         &
    !$omp redutcion(min:dmin,vmin)
    do k=1,Ngrid
       do j=1,Ngrid
          do i=1,Ngrid
             d = density(i, j, k)
             v = velocity(:, i, j, k)

             davg = davg + d
             dstd = dstd + d**2
             dmax = max(d, dmax)
             dmin = min(d, dmin)

             vavg = vavg + v
             vstd = vstd + v**2
             vmax = max(v, vmax)
             vmin = min(v, vmin)
          enddo
       enddo
    enddo
    !$omp end parallel do

    ! Finish computing statistics
    davg = davg / ncell
    dstd = sqrt(dstd / ncell - davg**2)
    vavg = vavg / ncell
    vstd = sqrt(vstd / ncell - vavg**2)

    write(*,*) "redshift: ", z
    write(*,*) "density:  ", real((/ davg, dstd, dmin, dmax /))
    write(*,*) "vx:       ", real((/ vavg(1), vstd(1), vmin(1), vmax(1) /))
    write(*,*) "vy:       ", real((/ vavg(2), vstd(2), vmin(2), vmax(2) /))
    write(*,*) "vz:       ", real((/ vavg(3), vstd(3), vmin(3), vmax(3) /))

    return


  contains


    ! Mass assignment

    subroutine ngp_mass_assignment(indx)
      ! Subroutine arguments
      integer(4), dimension(2,3) :: indx
      ! Local variables
      integer(4) :: i,j,k
      integer(4) :: ii,jj,kk
      integer(8) :: ip
      real(8)    :: x,y,z

      ! Loop over cells in domain
      do k=indx(1, 3),indx(2, 3)
      do j=indx(1, 2),indx(2, 2)
      do i=indx(1, 1),indx(2, 1)
         ! Head of chain
         ip = hoc_dm(i, j, k)

         do while(ip > 0)
            ! Particle
            x   = mod(dm_ip(1, ip) * x_dm_gr - 0.5 + ngridr, ngridr)
            y   = mod(dm_ip(2, ip) * x_dm_gr - 0.5 + ngridr, ngridr)
            z   = mod(dm_ip(3, ip) * x_dm_gr - 0.5 + ngridr, ngridr)

            ! Assign to grid
            ii = 1 + int(x)
            jj = 1 + int(y)
            kk = 1 + int(z)

            density(ii, jj, kk) = density(ii, jj, kk) + mass

            ! Next particle
            ip = ll_dm(ip)
         enddo
      enddo
      enddo
      enddo

      return
    end subroutine ngp_mass_assignment

    subroutine cic_mass_assignment(indx)
      ! Subroutine arguments
      integer(4), dimension(2,3) :: indx
      ! Local variables
      integer(4) :: i,j,k
      integer(4) :: i1,i2,j1,j2,k1,k2
      integer(8) :: ip
      real(8)    :: x,y,z
      real(8)    :: dx1,dx2,dy1,dy2,dz1,dz2

      ! Loop over cells in domain
      do k=indx(1, 3),indx(2, 3)
      do j=indx(1, 2),indx(2, 2)
      do i=indx(1, 1),indx(2, 1)
         ! Head of chain
         ip = hoc_dm(i, j, k)

         do while (ip > 0)
            ! Particle
            x   = mod(dm_ip(1, ip) * x_dm_gr - 0.5 + ngridr, ngridr)
            y   = mod(dm_ip(2, ip) * x_dm_gr - 0.5 + ngridr, ngridr)
            z   = mod(dm_ip(3, ip) * x_dm_gr - 0.5 + ngridr, ngridr)

            ! Assign CIC weights
            i1  = 1 + int(x)
            i2  = 1 + mod(i1, Ngrid)
            dx1 = i1 - x
            dx2 = 1 - dx1
            j1  = 1 + int(y)
            j2  = 1 + mod(j1, Ngrid)
            dy1 = j1 - y
            dy2 = 1 - dy1
            k1  = 1 + int(z)
            k2  = 1 + mod(k1, Ngrid)
            dz1 = k1 - z
            dz2 = 1 - dz1

            ! Assign mass to dm density field
            density(i1, j1, k1) = density(i1, j1, k1) + mass * dx1 * dy1 * dz1
            density(i2, j1, k1) = density(i2, j1, k1) + mass * dx2 * dy1 * dz1
            density(i1, j2, k1) = density(i1, j2, k1) + mass * dx1 * dy2 * dz1
            density(i2, j2, k1) = density(i2, j2, k1) + mass * dx2 * dy2 * dz1
            density(i1, j1, k2) = density(i1, j1, k2) + mass * dx1 * dy1 * dz2
            density(i2, j1, k2) = density(i2, j1, k2) + mass * dx2 * dy1 * dz2
            density(i1, j2, k2) = density(i1, j2, k2) + mass * dx1 * dy2 * dz2
            density(i2, j2, k2) = density(i2, j2, k2) + mass * dx2 * dy2 * dz2

            ! Next particle
            ip = ll_dm(ip)
         enddo
      enddo
      enddo
      enddo

      return
    end subroutine cic_mass_assignment

    subroutine tsc_mass_assignment(indx)
      ! Subroutine arguments
      integer(4), dimension(2,3) :: indx
      ! Local variables
      integer(4) :: i,j,k
      integer(4) :: a,b,c
      integer(8) :: ip
      real(8)    :: x,y,z,dx,dy,dz
      integer(4), dimension(-1:1) :: ii,jj,kk
      real(8),    dimension(-1:1) :: wx,wy,wz

      ! Loop over cells in domain
      do k=indx(1, 3),indx(2, 3)
      do j=indx(1, 2),indx(2, 2)
      do i=indx(1, 1),indx(2, 1)
         ! Head of chain
         ip = hoc_dm(i, j, k)

         do while (ip > 0)
            ! Particle
            x   = mod(dm_ip(1, ip) * x_dm_gr - 0.5 + ngridr, ngridr)
            y   = mod(dm_ip(2, ip) * x_dm_gr - 0.5 + ngridr, ngridr)
            z   = mod(dm_ip(3, ip) * x_dm_gr - 0.5 + ngridr, ngridr)

            ! Indices
            ii( 0) = 1 + int(x)
            ii(-1) = 1 + mod(ii(0) - 2 + Ngrid, Ngrid)
            ii( 1) = 1 + mod(ii(0)            , Ngrid)
            jj( 0) = 1 + int(y)
            jj(-1) = 1 + mod(jj(0) - 2 + Ngrid, Ngrid)
            jj( 1) = 1 + mod(jj(0)            , Ngrid)
            kk( 0) = 1 + int(z)
            kk(-1) = 1 + mod(kk(0) - 2 + Ngrid, Ngrid)
            kk( 1) = 1 + mod(kk(0)            , Ngrid)

            ! Weights
            dx     = x - (ii(0) - 0.5)
            wx( 0) = 0.75 - dx**2
            wx(-1) = 0.5*(1.5 - abs(dx + 1))**2
            wx( 1) = 0.5*(1.5 - abs(dx - 1))**2
            dy     = y - (jj(0) - 0.5)
            wy( 0) = 0.75 - dy**2
            wy(-1) = 0.5*(1.5 - abs(dy + 1))**2
            wy( 1) = 0.5*(1.5 - abs(dy - 1))**2
            dz     = z - (kk(0) - 0.5)
            wz( 0) = 0.75 - dz**2
            wz(-1) = 0.5*(1.5 - abs(dz + 1))**2
            wz( 1) = 0.5*(1.5 - abs(dz - 1))**2

            ! Add mass to density field
            do c=-1,1
            do b=-1,1
            do a=-1,1
               density(ii(a), jj(b), kk(c)) = density(ii(a), jj(b), kk(c)) &
                                            + mass * wx(a) * wy(b) * wz(c)
            enddo
            enddo
            enddo

            ! Next particle
            ip = ll_dm(ip)
         enddo
      enddo
      enddo
      enddo

      return
    end subroutine tsc_mass_assignment

    ! Velocity assignment

    subroutine ngp_vel_assignment(indx)
      ! Subroutine arguments
      integer(4), dimension(2,3) :: indx
      ! Local variables
      integer(4) :: i,j,k
      integer(4) :: ii,jj,kk
      integer(4) :: i1,i2,i3,idx(3)
      integer(8) :: ip
      real(8)    :: x,y,z,v(3)

      ! Loop over cells in domain
      do k=indx(1, 3),indx(2, 3)
      do j=indx(1, 2),indx(2, 2)
      do i=indx(1, 1),indx(2, 1)
         ! Head of chain
         ip = hoc_dm(i, j, k)

         do while (ip > 0)
            ! Particle
            x = mod(dm_ip(1, ip) * x_dm_gr - 0.5 + ngridr, ngridr)
            y = mod(dm_ip(2, ip) * x_dm_gr - 0.5 + ngridr, ngridr)
            z = mod(dm_ip(3, ip) * x_dm_gr - 0.5 + ngridr, ngridr)
            v = dm_ip(4:6, ip) * vel_unit / 1D5

            ! Assign velocity to grid
            ii = 1 + int(x)
            jj = 1 + int(y)
            kk = 1 + int(z)

            if (density(ii, jj, kk) > 0) then
               velocity(:, ii, jj, kk) = velocity(:, ii, jj, kk) &
                    + v * mass / density(ii,jj,kk)
            endif

            ! Next particle
            ip = ll_dm(ip)
         enddo
      enddo
      enddo
      enddo

      return
    end subroutine ngp_vel_assignment

    subroutine cic_vel_assignment(indx)
      ! Subroutine arguments
      integer(4), dimension(2,3) :: indx
      ! Local variables
      integer(4) :: i,j,k
      integer(4) :: i1,i2,j1,j2,k1,k2
      integer(8) :: ip
      real(8)    :: x,y,z,v(3)
      real(8)    :: dx1,dx2,dy1,dy2,dz1,dz2

      ! Loop over cells in domain
      do k=indx(1, 3),indx(2, 3)
      do j=indx(1, 2),indx(2, 2)
      do i=indx(1, 1),indx(2, 1)
         ! Head of chain
         ip = hoc_dm(i, j, k)

         do while (ip > 0)
            ! Particle
            x = mod(dm_ip(1, ip) * x_dm_gr - 0.5 + ngridr, ngridr)
            y = mod(dm_ip(2, ip) * x_dm_gr - 0.5 + ngridr, ngridr)
            z = mod(dm_ip(3, ip) * x_dm_gr - 0.5 + ngridr, ngridr)
            v = dm_ip(4:6, ip) * vel_unit / 1D5

            ! Assign CIC weights
            i1  = 1 + int(x)
            i2  = 1 + mod(i1, Ngrid)
            dx1 = i1 - x
            dx2 = 1 - dx1
            j1  = 1 + int(y)
            j2  = 1 + mod(j1, Ngrid)
            dy1 = j1 - y
            dy2 = 1 - dy1
            k1  = 1 + int(z)
            k2  = 1 + mod(k1, Ngrid)
            dz1 = k1 - z
            dz2 = 1 - dz1

            ! Assign velocity to field
            if (density(i1, j1, k1) > 0) then
               velocity(:, i1, j1, k1) = velocity(:, i1, j1, k1) &
                    + v * mass * dx1 * dy1 * dz1 / density(i1, j1, k1)
            endif
            if (density(i2, j1, k1) > 0) then
               velocity(:, i2, j1, k1) = velocity(:, i2, j1, k1) &
                    + v * mass * dx2 * dy1 * dz1 / density(i2, j1, k1)
            endif
            if (density(i1, j2, k1) > 0) then
               velocity(:, i1, j2, k1) = velocity(:, i1, j2, k1) &
                    + v * mass * dx1 * dy2 * dz1 / density(i1, j2, k1)
            endif
            if (density(i2, j2, k1) > 0) then
               velocity(:, i2, j2, k1) = velocity(:, i2, j2, k1) &
                    + v * mass * dx2 * dy2 * dz1 / density(i2, j2, k1)
            endif
            if (density(i1, j1, k2) > 0) then
               velocity(:, i1, j1, k2) = velocity(:, i1, j1, k2) &
                    + v * mass * dx1 * dy1 * dz2 / density(i1, j1, k2)
            endif
            if (density(i2, j1, k2) > 0) then
               velocity(:, i2, j1, k2) = velocity(:, i2, j1, k2) &
                    + v * mass * dx2 * dy1 * dz2 / density(i2, j1, k2)
            endif
            if (density(i1, j2, k2) > 0) then
               velocity(:, i1, j2, k2) = velocity(:, i1, j2, k2) &
                    + v * mass * dx1 * dy2 * dz2 / density(i1, j2, k2)
            endif
            if (density(i2, j2, k2) > 0) then
               velocity(:, i2, j2, k2) = velocity(:, i2, j2, k2) &
                    + v * mass * dx2 * dy2 * dz2 / density(i2, j2, k2)
            endif

            ! Next particle
            ip = ll_dm(ip)
         enddo
      enddo
      enddo
      enddo

      return
    end subroutine cic_vel_assignment

    subroutine tsc_vel_assignment(indx)
      ! Subroutine arguments
      integer(4), dimension(2,3) :: indx
      ! Local variables
      integer(4) :: i,j,k
      integer(4) :: a,b,c
      integer(8) :: ip
      real(8)    :: x,y,z,dx,dy,dz,v(3)
      integer(4), dimension(-1:1) :: ii,jj,kk
      real(8),    dimension(-1:1) :: wx,wy,wz

      ! Loop over cells in domain
      do k=indx(1, 3),indx(2, 3)
      do j=indx(1, 2),indx(2, 2)
      do i=indx(1, 1),indx(2, 1)
         ! Head of chain
         ip = hoc_dm(i, j, k)

         do while (ip > 0)
            ! Particle
            x   = mod(dm_ip(1, ip) * x_dm_gr - 0.5 + ngridr, ngridr)
            y   = mod(dm_ip(2, ip) * x_dm_gr - 0.5 + ngridr, ngridr)
            z   = mod(dm_ip(3, ip) * x_dm_gr - 0.5 + ngridr, ngridr)
            v   = dm_ip(4:6, ip) * vel_unit / 1D5

            ! Indices
            ii( 0) = 1 + int(x)
            ii(-1) = 1 + mod(ii(0) - 2 + Ngrid, Ngrid)
            ii( 1) = 1 + mod(ii(0)            , Ngrid)
            jj( 0) = 1 + int(y)
            jj(-1) = 1 + mod(jj(0) - 2 + Ngrid, Ngrid)
            jj( 1) = 1 + mod(jj(0)            , Ngrid)
            kk( 0) = 1 + int(z)
            kk(-1) = 1 + mod(kk(0) - 2 + Ngrid, Ngrid)
            kk( 1) = 1 + mod(kk(0)            , Ngrid)

            ! Weights
            dx     = x - (ii(0) - 0.5)
            wx( 0) = 0.75 - dx**2
            wx(-1) = 0.5*(1.5 - abs(dx + 1))**2
            wx( 1) = 0.5*(1.5 - abs(dx - 1))**2
            dy     = y - (jj(0) - 0.5)
            wy( 0) = 0.75 - dy**2
            wy(-1) = 0.5*(1.5 - abs(dy + 1))**2
            wy( 1) = 0.5*(1.5 - abs(dy - 1))**2
            dz     = z - (kk(0) - 0.5)
            wz( 0) = 0.75 - dz**2
            wz(-1) = 0.5*(1.5 - abs(dz + 1))**2
            wz( 1) = 0.5*(1.5 - abs(dz - 1))**2

            ! Add mass to density field
            do c=-1,1
            do b=-1,1
            do a=-1,1
               if (density(ii(a), jj(b), kk(c)) > 0) then
                  velocity(:, ii(a), jj(b), kk(c)) =      &
                       velocity(:, ii(a), jj(b), kk(c))   &
                       + v * mass * wx(a) * wy(b) * wz(c) &
                       / density(ii(a), jj(b), kk(c))
               endif
            enddo
            enddo
            enddo

            ! Next particle
            ip = ll_dm(ip)
         enddo
      enddo
      enddo
      enddo

      return
    end subroutine tsc_vel_assignment


  end subroutine calc_density_velocity_fields


  subroutine calc_ion_t21_fields(z, density, zreion, ion, t21)
    ! Subroutine arguments
    real(c_double), intent(in)    :: z
    real(c_double), intent(in)    :: density(:,:,:),zreion(:,:,:)
    real(c_double), intent(inout) :: ion(:,:,:),t21(:,:,:)
    ! Local variables
    integer(4) :: i,j,k
    real(8)    :: ncell,t0,h0,ob,om
    real(8)    :: x,xavg,xstd,xmin,xmax
    real(8)    :: t,tavg,tstd,tmin,tmax

    ! Unpack parameters
    h0    = hubble0
    ob    = omegab
    om    = omegam
    ncell = real(Ngrid, kind=8)**3

    ! Calculate global temperature
    t0 = T0_of_z(z, h0, ob, om)

    ! Compute ionization field, 21cm field and statistics
    xavg = 0
    xstd = 0
    xmin =  huge(xmin)
    xmax = -huge(xmax)
    tavg = 0
    tstd = 0
    tmin =  huge(tmin)
    tmax = -huge(tmax)
    !$omp parallel do                      &
    !$omp default(shared)                  &
    !$omp private(i,j,k,x,t)               &
    !$omp reduction(+:xavg,xstd,tavg,tstd) &
    !$omp reduction(max:xmax,tmax)         &
    !$omp reduction(min:xmin,tmin)
    do k=1,Ngrid
       do j=1,Ngrid
          do i=1,Ngrid
             if (zreion(i, j, k) > z) then
                x = 1
                t = 0
             else
                x = 0
                t = t0 * density(i, j, k)
             endif

             ion(i, j, k) = x
             t21(i, j, k) = t

             xavg = xavg + x
             xstd = xstd + x**2
             xmax = max(x, xmax)
             xmin = min(x, xmin)

             tavg = tavg + t
             tstd = tstd + t**2
             tmax = max(t, tmax)
             tmin = min(t, tmin)
          enddo
       enddo
    enddo
    !$omp end parallel do

    ! Write out statistics
    xavg = xavg / ncell
    xstd = sqrt(xstd / ncell - xavg**2)
    tavg = tavg / ncell
    tstd = sqrt(tstd / ncell - tavg**2)
    write(*,*) "ion: ", real((/ xavg, xstd, xmin, xmax /))
    write(*,*) "T21: ", real((/ tavg, tstd, tmin, tmax /))

    return


  contains


    pure function T0_of_z(z, h0, ob, om)
      ! Function arguments
      real(8), intent(in) :: z,h0,ob,om
      real(8)             :: T0_of_z

      ! Unpack cosmology
      T0_of_z = 38.6D0 * h0 * (ob / 0.045D0) &
           * sqrt(0.27D0 / om * (1 + z) / 10D0)
      return
    end function T0_of_z


  end subroutine calc_ion_t21_fields


end module field_tools
