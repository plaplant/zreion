! MKL
include 'mkl_vsl.f90'

module ic_tools
  ! OpenMP
  use omp_lib

  ! MKL
  use mkl_vsl_type
  use mkl_vsl

  ! Imports for C linking
  use iso_c_binding

  ! general tools
  use general_tools

  ! FFT tools
  use fft_tools

  ! Cosmo tools
  use cosmo_tools

  ! Default
  implicit none


contains


  subroutine generate_ics(field, pk_lin, seed)
    ! Subroutine arguments
    real(c_double), intent(inout) :: field(:,:,:)
    real(c_double), intent(in)    :: pk_lin(:,:)
    integer(c_int), intent(in)    :: seed
    ! Local variables
    integer(4) :: i,j,k,Nx,Ny,Nz
    integer(4) :: brng,method,errcode
    real(8)    :: a,b,sigma,ngrid
    real(8)    :: d,davg,dstd,dmax,dmin
    real(8)    :: kx,ky,kz,kr,pk,x
    type(vsl_stream_state)  :: stream
    integer(4), allocatable :: iseed(:)
    real(8),    allocatable :: xseed(:)

    ! Unpack arguments
    Nx    = size(field, dim=1) - 2  ! FFT padding
    Ny    = size(field, dim=2)
    Nz    = size(field, dim=3)
    ngrid = int(Nx, kind=8) * int(Ny, kind=8) * int(Nz, kind=8)

    ! Allocate
    allocate(iseed(Nz))
    allocate(xseed(Nz))

    ! Generate GRF
    iseed(1) = seed
    a        = 0
    b        = 1
    brng     = VSL_BRNG_MT19937
    method   = VSL_RNG_METHOD_UNIFORM_STD_ACCURATE
    errcode  = vslnewstream(stream, brng, iseed(1))
    errcode  = vdrnguniform(method, stream, Nz, xseed, a, b)
    iseed    = int(xseed * huge(0))

    ! Generate Gaussian random numbers
    brng   = VSL_BRNG_MT19937
    method = VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2
    a      = 0
    sigma  = 1

    !$omp parallel do                 &
    !$omp default(shared)             &
    !$omp private(j,k,stream,errcode)
    do k=1,Nz
       ! Initialize stream
       errcode = vslnewstream(stream, brng, iseed(k))

       ! Generate random numbers
       do j=1,Ny
          errcode = vdrnggaussian(method, stream, Nx, field(1:Nx, j, k), a, &
               sigma)
       enddo

       ! Delete stream
       errcode = vsldeletestream(stream)
    enddo
    !$omp end parallel do

    ! Check Gaussian random field
    davg = 0
    dstd = 0
    dmax = -huge(dmax)
    dmin =  huge(dmin)

    !$omp parallel do            &
    !$omp default(shared)        &
    !$omp private(i,j,k,d)       &
    !$omp reduction(+:davg,dstd) &
    !$omp reduction(max:dmax)    &
    !$omp reduction(min:dmin)
    do k=1,Nz
       do j=1,Ny
          do i=1,Nx
             d    = field(i,j,k)
             davg = davg + d
             dstd = dstd + d**2
             dmax = max(d,dmax)
             dmin = min(d,dmin)
          enddo
       enddo
    enddo
    !$omp end parallel do

    ! Finish computing mean and std
    davg  = davg / ngrid
    dstd  = sqrt(dstd / ngrid - davg**2)
    write(*,*) "GRF: ",real((/ davg, dstd, dmin, dmax /))

    ! Compute FFT
    call fft_3d(field)

    ! Interpolate from linear P(k)
    !$omp parallel do          &
    !$omp default(shared)      &
    !$omp private(i,j,k,pk)    &
    !$omp private(kx,ky,kz,kr)
    do k=1,Nz
       if (k <= Nz / 2 + 1) then
          kz = 2 * pi / box * (k - 1)
       else
          kz = 2 * pi / box * (k - 1 - Nz)
       endif

       do j=1,Ny
          if (j <= Ny / 2 + 1) then
             ky = 2 * pi / box * (j - 1)
          else
             ky = 2 * pi / box * (j - 1 - Ny)
          endif

          do i=1,Nx+2,2
             kx = 2 * pi / box * ((i - 1) / 2)
             kr = sqrt(kx**2 + ky**2 + kz**2)

             if (kr > 0) then
                ! Interpolate from linear theory
                pk = interpolate(kr, pk_lin, 1, 2, "log")
                ! Save in perturbation field
                field(i:i+1, j, k) = sqrt(pk) * field(i:i+1, j, k)
             else
                field(i:i+1, j, k) = 0
             endif
          enddo
       enddo
    enddo
    !$omp end parallel do

    return
  end subroutine generate_ics

end module ic_tools
