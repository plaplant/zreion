module zreion_tools
  ! OpenMP
  use omp_lib

  ! Imports for C linking
  use iso_c_binding

  ! FFT tools
  use fft_tools

  ! Cosmology
  use cosmo_tools

  ! General tools
  use general_tools

  ! Default
  implicit none


contains


  subroutine calc_zreion(density, zreion)
    ! Subroutine arguments
    real(c_double), intent(in)    :: density(:,:,:)
    real(c_double), intent(inout) :: zreion(:,:,:)
    ! Local variables
    integer(4) :: i,j,k
    real(8)    :: kx,ky,kz,kr,Ak
    real(8)    :: wcell,wsmooth,x
    real(8)    :: z,zavg,zstd,zmax,zmin

    ! Unpack parameters
    Ak = 2 * pi / real(Ngrid, kind=8)

    ! Copy density field to zreion
    !$omp parallel do     &
    !$omp default(shared) &
    !$omp private(i,j,k)
    do k=1,Ngrid
       do j=1,Ngrid
          do i=1,Ngrid
             zreion(i,j,k) = density(i,j,k) - 1
          enddo
       enddo
    enddo
    !$omp end parallel do

    ! Forward FFT of density field
    call fft_3d(zreion)

    ! Apply bias factor in Fourier space
    !$omp parallel do                &
    !$omp default(shared)            &
    !$omp private(i,j,k,kx,ky,kz,kr) &
    !$omp private(x,wcell,wsmooth)
    do k=1,Ngrid
       if (k <= Ngrid / 2 + 1) then
          kz = Ak * (k - 1)
       else
          kz = Ak * (k - 1 - Ngrid)
       endif

       do j=1,Ngrid
          if (j <= Ngrid / 2 + 1) then
             ky = Ak * (j - 1)
          else
             ky = Ak * (j - 1 - Ngrid)
          endif

          do i=1,Ngrid+2,2
             kx = Ak * ((i - 1) / 2)
             kr = sqrt(kx**2 + ky**2 + kz**2)

             ! Deconvolve cell
             wcell = (sinc(kx / 2) * sinc(ky / 2) * sinc(kz / 2))**2

             ! Smooth with top hat window
             x = kr * Rsmooth_zre
             wsmooth = tophat(x)

             ! Apply bias relation
             zreion(i:i+1, j, k) = zreion(i:i+1, j, k) * &
                  bias(kr, b0_zre, kb_zre, alpha_zre) * wsmooth / wcell
          enddo
       enddo
    enddo
    !$omp end parallel do

    ! Inverse FFT of zreion
    call ifft_3d(zreion)

    ! Finish calculating field and statistics
    zavg = 0
    zstd = 0
    zmax = -huge(zmax)
    zmin =  huge(zmin)

    !$omp parallel do            &
    !$omp default(shared)        &
    !$omp private(i,j,k,z)       &
    !$omp reduction(+:zavg,zstd) &
    !$omp reduction(max:zmax)    &
    !$omp reduction(min:zmin)
    do k=1,Ngrid
       do j=1,Ngrid
          do i=1,Ngrid
             z = zmean_zre + (1 + zmean_zre) * zreion(i, j, k)
             zreion(i, j, k) = z

             zavg = zavg + z
             zstd = zstd + z**2
             zmax = max(z, zmax)
             zmin = min(z, zmin)
          enddo
       enddo
    enddo
    !$omp end parallel do

    ! Write out zreion properties
    x = real(Ngrid, kind=8)**3
    zavg = zavg / x
    zstd = sqrt(zstd / x - zavg**2)
    write(*,*) "zreion: ", real((/ zavg, zstd, zmin, zmax /))

    return

  contains

    pure function bias(kr, b0, kb, alpha)
      ! Function parameters
      real(8), intent(in) :: kr,b0,kb,alpha
      real(8)             :: bias

      bias = b0 / (1 + kr / kb)**alpha
      return
    end function bias

  end subroutine calc_zreion


end module zreion_tools
