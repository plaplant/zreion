module particle_tools
  ! OpenMP
  use omp_lib

  ! Imports for C linking
  use iso_c_binding

  ! Cosmology tools
  use cosmo_tools

  ! FFT tools
  use fft_tools

  ! Default
  implicit none

  ! Particle arrays
  real(c_float),   allocatable, dimension(:,:,:,:), target  :: dm
  real(c_float),                dimension(:,:),     pointer :: dm_ip
  integer(c_long), allocatable, dimension(:)                :: ll_dm
  integer(c_long), allocatable, dimension(:,:)              :: dm_cpu
  integer(c_long), allocatable, dimension(:,:,:)            :: hoc_dm,hoc_cpu


contains


  subroutine init_dm
    ! Local variables
    integer(4) :: icpu,i,j,k
    integer(8) :: Npart,npt

    ! Optionally deallocate arrays
    if (allocated(dm)) then
       deallocate(dm)
       deallocate(ll_dm)
       deallocate(hoc_dm)
       deallocate(hoc_cpu)
    endif

    ! Allocate arrays
    Npart = int(Ndm, kind=8)**3
    allocate(dm(6, Ndm, Ndm, Ndm))
    allocate(dm_cpu(2, Ncpu))
    allocate(ll_dm(Npart))
    allocate(hoc_dm(Ngrid, Ngrid, Ngrid))
    allocate(hoc_cpu(Ngrid, Ngrid, Ncpu))

    ! Initialize in parallel
    npt = Npart / Ncpu + min(1, mod(Npart, Ncpu))
    !$omp parallel        &
    !$omp default(shared) &
    !$omp private(i,j,k)
    !$omp do
    do k=1,Ndm
       do j=1,Ndm
          do i=1,Ndm
             dm(:,i,j,k) = 0
          enddo
       enddo
    enddo
    !$omp end do
    !$omp do
    do icpu=1,Ncpu
       dm_cpu(1,icpu) = 1 + (icpu - 1) * npt
       dm_cpu(2,icpu) = min(icpu * npt, Npart)
    enddo
    !$omp end do
    !$omp do
    do icpu=1,Ncpu
       ll_dm(dm_cpu(1,icpu):dm_cpu(2,icpu)) = 0
    enddo
    !$omp end do
    !$omp do
    do k=1,Ngrid
       hoc_dm(:,:,k) = 0
    enddo
    !$omp end do
    !$omp do
    do k=1,Ncpu
       hoc_cpu(:,:,k) = 0
    enddo
    !$omp end do
    !$omp end parallel

    ! Associate pointer
    dm_ip(1:6, 1:Npart) => dm

    return
  end subroutine init_dm


  subroutine calc_2lpt(z, delta1)
    ! Subroutine arguments
    real(c_double), intent(in) :: z
    real(c_double), intent(in) :: delta1(:,:,:)
    ! Local variables
    integer(4) :: i,j,k
    real(8)    :: kx,ky,kz,ksq
    real(8)    :: wx,wy,wz
    real(8), allocatable  :: gradphi(:,:,:,:),delta2(:,:,:)

    ! Allocate and initialize in parallel
    allocate(gradphi(Ndm + 2, Ndm, Ndm, 3))
    allocate(delta2(Ndm + 2, Ndm, Ndm))

    !$omp parallel        &
    !$omp default(shared) &
    !$omp private(k)
    !$omp do
    do k=1,Ndm
       gradphi(:, :, k, 1) = 0
    enddo
    !$omp end do
    !$omp do
    do k=1,Ndm
       gradphi(:, :, k, 2) = 0
    enddo
    !$omp end do
    !$omp do
    do k=1,Ndm
       gradphi(:, :, k, 3) = 0
    enddo
    !$omp end do
    !$omp do
    do k=1,Ndm
       delta2(:, :, k) = 0
    enddo
    !$omp end do
    !$omp end parallel

    ! Calcualte growth factors
    call calc_cosmo_vars(z)

    ! Calcualte Grad phi_1 fields
    !$omp parallel do           &
    !$omp default(shared)       &
    !$omp private(i,j,k)        &
    !$omp private(kx,ky,kz,ksq) &
    !$omp private(wx,wy,wz)
    do k=1,Ndm
       if (k <= Ndm / 2 + 1) then
          kz = 2 * pi / real(Ndm, kind=8) * (k - 1)
       else
          kz = 2 * pi / real(Ndm, kind=8) * (k - 1 - Ndm)
       endif

       do j=1,Ndm
          if (j <= Ndm / 2 + 1) then
             ky = 2 * pi / real(Ndm, kind=8) * (j - 1)
          else
             ky = 2 * pi / real(Ndm, kind=8) * (j - 1 - Ndm)
          endif

          do i=1,Ndm+2,2
             kx = 2 * pi / real(Ndm, kind=8) * ((i - 1) / 2)

             ksq = kx**2 + ky**2 + kz**2

             ! Force kernels
             if (ksq > 0) then
                wx = -kx / ksq
                wy = -ky / ksq
                wz = -kz / ksq
             else
                wx = 0
                wy = 0
                wz = 0
             endif

             ! Complex multiply
             gradphi(i    , j, k, 1) = -wx * delta1(i + 1, j, k)
             gradphi(i + 1, j, k, 1) =  wx * delta1(i    , j, k)
             gradphi(i    , j, k, 2) = -wy * delta1(i + 1, j, k)
             gradphi(i + 1, j, k, 2) =  wy * delta1(i    , j, k)
             gradphi(i    , j, k, 3) = -wz * delta1(i + 1, j, k)
             gradphi(i + 1, j, k, 3) =  wz * delta1(i    , j, k)
          enddo
       enddo
    enddo
    !$omp end parallel do

    ! Inverse FFT
    call ifft_3d(gradphi(:, :, :, 1))
    call ifft_3d(gradphi(:, :, :, 2))
    call ifft_3d(gradphi(:, :, :, 3))

    ! Save in particle array
    !$omp parallel do     &
    !$omp default(shared) &
    !$omp private(i,j,k)
    do k=1,Ndm
       do j=1,Ndm
          do i=1,Ndm
             dm(1, i, j, k) = -Dfactor(1) * gradphi(i, j, k, 1)
             dm(2, i, j, k) = -Dfactor(1) * gradphi(i, j, k, 2)
             dm(3, i, j, k) = -Dfactor(1) * gradphi(i, j, k, 3)
             dm(4, i, j, k) = -vfactor(1) * gradphi(i, j, k, 1)
             dm(5, i, j, k) = -vfactor(1) * gradphi(i, j, k, 2)
             dm(6, i, j, k) = -vfactor(1) * gradphi(i, j, k, 3)
          enddo
       enddo
    enddo
    !$omp end parallel do

    ! Calculate Grad^2 phi_ii fields
    !$omp parallel do           &
    !$omp default(shared)       &
    !$omp private(i,j,k)        &
    !$omp private(kx,ky,kz,ksq) &
    !$omp private(wx,wy,wz)
    do k=1,Ndm
       if (k <= Ndm / 2 + 1) then
          kz = 2 * pi / real(Ndm, kind=8) * (k - 1)
       else
          kz = 2 * pi / real(Ndm, kind=8) * (k - 1 - Ndm)
       endif

       do j=1,Ndm
          if (j <= Ndm / 2 + 1) then
             ky = 2 * pi / real(Ndm, kind=8) * (j - 1)
          else
             ky = 2 * pi / real(Ndm, kind=8) * (j - 1 - Ndm)
          endif

          do i=1,Ndm+2,2
             kx = 2 * pi / real(Ndm, kind=8) * ((i - 1) / 2)

             ksq = kx**2 + ky**2 + kz**2

             ! Calc kernels
             if (ksq > 0) then
                wx = -kx**2 / ksq
                wy = -ky**2 / ksq
                wz = -kz**2 / ksq
             else
                wx = 0
                wy = 0
                wz = 0
             endif

             ! Complex multiply
             gradphi(i:i+1, j, k, 1) = wx * delta1(i:i+1, j, k)
             gradphi(i:i+1, j, k, 2) = wy * delta1(i:i+1, j, k)
             gradphi(i:i+1, j, k, 3) = wz * delta1(i:i+1, j, k)
          enddo
       enddo
    enddo
    !$omp end parallel do

    ! Inverse FFT Grad^2 phi_ii fields
    call ifft_3d(gradphi(:, :, :, 1))
    call ifft_3d(gradphi(:, :, :, 2))
    call ifft_3d(gradphi(:, :, :, 3))

    ! Calculate delta2
    !$omp parallel do     &
    !$omp default(shared) &
    !$omp private(i,j,k)
    do k=1,Ndm
       do j=1,Ndm
          do i=1,Ndm
             delta2(i, j, k) = gradphi(i, j, k, 1) * gradphi(i, j, k, 2) &
                             + gradphi(i, j, k, 2) * gradphi(i, j, k, 3) &
                             + gradphi(i, j, k, 3) * gradphi(i, j, k, 1)
          enddo
       enddo
    enddo
    !$omp end parallel do

    ! Calculate Grad^2 phi_ij fields
    !$omp parallel do           &
    !$omp default(shared)       &
    !$omp private(i,j,k)        &
    !$omp private(kx,ky,kz,ksq) &
    !$omp private(wx,wy,wz)
    do k=1,Ndm
       if (k <= Ndm / 2 + 1) then
          kz = 2 * pi / real(Ndm, kind=8) * (k - 1)
       else
          kz = 2 * pi / real(Ndm, kind=8) * (k - 1 - Ndm)
       endif

       do j=1,Ndm
          if (j <= Ndm / 2 + 1) then
             ky = 2 * pi / real(Ndm, kind=8) * (j - 1)
          else
             ky = 2 * pi / real(Ndm, kind=8) * (j - 1 - Ndm)
          endif

          do i=1,Ndm+2,2
             kx = 2 * pi / real(Ndm, kind=8) * ((i - 1) / 2)

             ksq = kx**2 + ky**2 + kz**2

             ! Calc kernels
             if (ksq > 0) then
                wx = -kx * ky / ksq
                wy = -ky * kz / ksq
                wz = -kz * kx / ksq
             else
                wx = 0
                wy = 0
                wz = 0
             endif

             ! Complex multiply
             gradphi(i:i+1, j, k, 1) = wx * delta1(i:i+1, j, k)
             gradphi(i:i+1, j, k, 2) = wy * delta1(i:i+1, j, k)
             gradphi(i:i+1, j, k, 3) = wz * delta1(i:i+1, j, k)
          enddo
       enddo
    enddo
    !$omp end parallel do

    ! Inverse FFT Grad^2 phi_ij fields
    call ifft_3d(gradphi(:, :, :, 1))
    call ifft_3d(gradphi(:, :, :, 2))
    call ifft_3d(gradphi(:, :, :, 3))

    ! Add to delta2
    !$omp parallel do     &
    !$omp default(shared) &
    !$omp private(i,j,k)
    do k=1,Ndm
       do j=1,Ndm
          do i=1,Ndm
             delta2(i, j, k) = delta2(i, j, k) &
                             - gradphi(i, j, k, 1)**2 &
                             - gradphi(i, j, k, 2)**2 &
                             - gradphi(i, j, k, 3)**2
          enddo
       enddo
    enddo
    !$omp end parallel do

    ! Forward FFT delta2 field
    call fft_3d(delta2)

    ! Calculate Grad phi_2
    !$omp parallel do           &
    !$omp default(shared)       &
    !$omp private(i,j,k)        &
    !$omp private(kx,ky,kz,ksq) &
    !$omp private(wx,wy,wz)
    do k=1,Ndm
       if (k <= Ndm / 2 + 1) then
          kz = 2 * pi / real(Ndm, kind=8) * (k - 1)
       else
          kz = 2 * pi / real(Ndm, kind=8) * (k - 1 - Ndm)
       endif

       do j=1,Ndm
          if (j <= Ndm / 2 + 1) then
             ky = 2 * pi / real(Ndm, kind=8) * (j - 1)
          else
             ky = 2 * pi / real(Ndm, kind=8) * (j - 1 - Ndm)
          endif

          do i=1,Ndm+2,2
             kx = 2 * pi / real(Ndm, kind=8) * ((i - 1) / 2)

             ksq = kx**2 + ky**2 + kz**2

             ! Kernels
             if (ksq > 0) then
                wx = -kx / ksq
                wy = -ky / ksq
                wz = -kz / ksq
             else
                wx = 0
                wy = 0
                wz = 0
             endif

             ! Complex multiply
             gradphi(i  , j, k, 1) = -wx * delta2(i+1, j, k)
             gradphi(i+1, j, k, 1) =  wx * delta2(i  , j, k)
             gradphi(i  , j, k, 2) = -wy * delta2(i+1, j, k)
             gradphi(i+1, j, k, 2) =  wy * delta2(i  , j, k)
             gradphi(i  , j, k, 3) = -wz * delta2(i+1, j, k)
             gradphi(i+1, j, k, 3) =  wz * delta2(i  , j, k)
          enddo
       enddo
    enddo
    !$omp end parallel do

    ! Inverse FFT Grad phi_2 fields
    call ifft_3d(gradphi(:, :, :, 1))
    call ifft_3d(gradphi(:, :, :, 2))
    call ifft_3d(gradphi(:, :, :, 3))

    ! Save in particle array
    !$omp parallel do     &
    !$omp default(shared) &
    !$omp private(i,j,k)
    do k=1,Ndm
       do j=1,Ndm
          do i=1,Ndm
             dm(1, i, j, k) = dm(1, i, j, k) &
                              + Dfactor(2) * gradphi(i, j, k, 1)
             dm(2, i, j, k) = dm(2, i, j, k) &
                              + Dfactor(2) * gradphi(i, j, k, 2)
             dm(3, i, j, k) = dm(3, i, j, k) &
                              + Dfactor(2) * gradphi(i, j, k, 3)
             dm(4, i, j, k) = dm(4, i, j, k) &
                              + vfactor(2) * gradphi(i, j, k, 1)
             dm(5, i, j, k) = dm(5, i, j, k) &
                              + vfactor(2) * gradphi(i, j, k, 2)
             dm(6, i, j, k) = dm(6, i, j, k) &
                              + vfactor(2) * gradphi(i, j, k, 3)
          enddo
       enddo
    enddo
    !$omp end parallel do

    ! Deallocate
    deallocate(gradphi)
    deallocate(delta2)

    return
  end subroutine calc_2lpt


  subroutine calc_dm
    ! Local variables
    integer(4) :: i,j,k
    real(8)    :: x(3),v(3),ndmr

    ! Unpack parameters
    ndmr = real(Ndm, kind=8)

    ! Compute particle positions + displacement
    !$omp parallel do        &
    !$omp default(shared)    &
    !$omp private(i,j,k,x,v)
    do k=1,Ndm
       do j=1,Ndm
          do i=1,Ndm
             x = mod(dm(1:3, i, j, k) + (/ i, j, k /) - 0.5 + ndmr, ndmr)
             v = dm(4:6, i, j, k) * time_unit

             dm(1:3, i, j, k) = x
             dm(4:6, i, j, k) = v
          enddo
       enddo
    enddo
    !$omp end parallel do

    return
  end subroutine calc_dm


  subroutine link_dm
    ! Local variables
    integer(4) :: icpu,i,j,k
    integer(4) :: i1,i2,i3,idx(3)
    integer(8) :: ip,next
    real(8)    :: x_dm_gr

    ! Unpack arguments
    x_dm_gr = real(Ngrid, kind=8) / real(Ndm, kind=8)

    ! Construct chaining lists in parallel
    !$omp parallel             &
    !$omp default(shared)      &
    !$omp private(icpu,ip,idx) &
    !$omp private(i,j,k,next)  &
    !$omp private(i1,i2,i3)
    !$omp do
    do icpu=1,Ncpu
       ! Init
       hoc_cpu(:, :, icpu) = 0

       ! Loop over particles
       do ip=dm_cpu(1,icpu),dm_cpu(2,icpu)
          j                   = 1 + mod(int(dm_ip(2, ip) * x_dm_gr), Ngrid)
          k                   = 1 + mod(int(dm_ip(3, ip) * x_dm_gr), Ngrid)
          ll_dm(ip)           = hoc_cpu(j, k, icpu)
          hoc_cpu(j, k, icpu) = ip
       enddo
    enddo
    !$omp end do
    !$omp do
    do k=1,Ngrid
       do j=1,Ngrid
          ! Init
          hoc_dm(:, j, k) = 0

          ! Loop over cpus
          do icpu=Ncpu,1,-1
             ! Head of chain
             ip = hoc_cpu(j, k, icpu)

             do while (ip > 0)
                i               = 1 + mod(int(dm_ip(1, ip) * x_dm_gr), Ngrid)
                next            = ll_dm(ip)
                ll_dm(ip)       = hoc_dm(i, j, k)
                hoc_dm(i, j, k) = ip
                ip              = next
             enddo
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel

    return
  end subroutine link_dm


end module particle_tools
