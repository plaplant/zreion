! MKL
include 'mkl_dfti.f90'

module fft_tools
  ! FFT
  use mkl_dfti

  ! Imports for C linking
  use iso_c_binding

  ! Default
  implicit none

  ! Interfaces for different precisions
  interface fft_2d
     module procedure fft_2d_s, fft_2d_d
  end interface fft_2d

  interface ifft_2d
     module procedure ifft_2d_s, ifft_2d_d
  end interface ifft_2d

  interface fft_3d
     module procedure fft_3d_s, fft_3d_d
  end interface fft_3d

  interface ifft_3d
     module procedure ifft_3d_s, ifft_3d_d
  end interface ifft_3d


contains


  subroutine fft_3d_s(a)
    ! Subroutine arguments
    real(c_float), dimension(:, :, :), target, intent(inout) :: a
    ! Local variables
    integer(4)  :: status,n1,n2,n3
    integer(8)  :: Nfft
    type(c_ptr) :: fft_ptr
    integer(4), dimension(3)             :: Lfft
    integer(4), dimension(4)             :: strides_in,strides_out
    type(dfti_descriptor), pointer       :: desc
    real(c_float), pointer, dimension(:) :: a_fft

    ! Get array dimensions
    n1 = size(a, dim=1) - 2  ! FFT padding
    n2 = size(a, dim=2)
    n3 = size(a, dim=3)

    ! Create pointer
    Nfft    = int(n1 + 2, kind=8) * int(n2, kind=8) * int(n3, kind=8)
    fft_ptr = c_loc(a)
    call c_f_pointer(fft_ptr, a_fft, [Nfft])

    ! Initialize
    Lfft        = (/ n1, n2, n3 /)
    strides_in  = (/ 0, 1, n1 + 2, (n1 + 2) * n2 /)
    strides_out = (/ 0, 1, n1 / 2 + 1, (n1 / 2 + 1) * n2 /)

    ! Initialize descriptor
    status = DftiCreateDescriptor(desc, DFTI_SINGLE, DFTI_REAL, 3, Lfft)
    status = DftiSetValue(desc, DFTI_CONJUGATE_EVEN_STORAGE, &
         DFTI_COMPLEX_COMPLEX)
    status = DftiSetValue(desc, DFTI_INPUT_STRIDES, strides_in)
    status = DftiSetValue(desc, DFTI_OUTPUT_STRIDES, strides_out)
    status = DftiCommitDescriptor(desc)

    ! Compute transform
    status = DftiComputeForward(desc, a_fft)

    ! Free descriptor
    status = DftiFreeDescriptor(desc)

    return
  end subroutine fft_3d_s


  subroutine ifft_3d_s(a)
    ! Subroutine arguments
    real(c_float), dimension(:, :, :), target, intent(inout) :: a
    ! Local variables
    integer(4)  :: status,n1,n2,n3
    integer(8)  :: Nfft
    real(8)     :: bscale
    type(c_ptr) :: fft_ptr
    integer(4), dimension(3)             :: Lfft
    integer(4), dimension(4)             :: strides_in,strides_out
    type(dfti_descriptor), pointer       :: desc
    real(c_float), pointer, dimension(:) :: a_fft

    ! Get array dimensions
    n1 = size(a, dim=1) - 2  ! FFT padding
    n2 = size(a, dim=2)
    n3 = size(a, dim=3)

    ! Create pointer
    Nfft    = int(n1 + 2, kind=8) * int(n2, kind=8) * int(n3, kind=8)
    fft_ptr = c_loc(a)
    call c_f_pointer(fft_ptr, a_fft, [Nfft])

    ! Initialize
    Lfft        = (/ n1, n2, n3 /)
    strides_in  = (/ 0, 1, n1 / 2 + 1, (n1 / 2 + 1) * n2 /)
    strides_out = (/ 0, 1, n1 + 2, (n1 + 2) * n2 /)
    bscale      = 1D0 / (int(n1, kind=8) * int(n2, kind=8) * int(n3, kind=8))

    ! Initialize descriptor
    status = DftiCreateDescriptor(desc, DFTI_SINGLE, DFTI_REAL, 3, Lfft)
    status = DftiSetValue(desc, DFTI_CONJUGATE_EVEN_STORAGE, &
         DFTI_COMPLEX_COMPLEX)
    status = DftiSetValue(desc, DFTI_INPUT_STRIDES, strides_in)
    status = DftiSetValue(desc, DFTI_OUTPUT_STRIDES, strides_out)
    status = DftiSetValue(desc, DFTI_BACKWARD_SCALE, bscale)
    status = DftiCommitDescriptor(desc)

    ! Compute transform
    status = DftiComputeBackward(desc, a_fft)

    ! Free descriptor
    status = DftiFreeDescriptor(desc)

    return
  end subroutine ifft_3d_s


  subroutine fft_3d_d(a)
    ! Subroutine arguments
    real(c_double), dimension(:, :, :), target, intent(inout) :: a
    ! Local variables
    integer(4)  :: status,n1,n2,n3
    integer(8)  :: Nfft
    type(c_ptr) :: fft_ptr
    integer(4), dimension(3)              :: Lfft
    integer(4), dimension(4)              :: strides_in,strides_out
    type(dfti_descriptor), pointer        :: desc
    real(c_double), pointer, dimension(:) :: a_fft

    ! Get array dimensions
    n1 = size(a, dim=1) - 2  ! FFT padding
    n2 = size(a, dim=2)
    n3 = size(a, dim=3)

    ! Create pointer
    Nfft    = int(n1 + 2, kind=8) * int(n2, kind=8) * int(n3, kind=8)
    fft_ptr = c_loc(a)
    call c_f_pointer(fft_ptr, a_fft, [Nfft])

    ! Initialize
    Lfft        = (/ n1, n2, n3/)
    strides_in  = (/ 0, 1, n1 + 2, (n1 + 2) * n2 /)
    strides_out = (/ 0, 1, n1 / 2 + 1, (n1 / 2 + 1) * n2 /)

    ! Initialize descriptor
    status = DftiCreateDescriptor(desc, DFTI_DOUBLE, DFTI_REAL, 3, Lfft)
    status = DftiSetValue(desc, DFTI_CONJUGATE_EVEN_STORAGE, &
         DFTI_COMPLEX_COMPLEX)
    status = DftiSetValue(desc, DFTI_INPUT_STRIDES, strides_in)
    status = DftiSetValue(desc, DFTI_OUTPUT_STRIDES, strides_out)
    status = DftiCommitDescriptor(desc)

    ! Compute transform
    status = DftiComputeForward(desc, a_fft)

    ! Free descriptor
    status = DftiFreeDescriptor(desc)

    return
  end subroutine fft_3d_d


  subroutine ifft_3d_d(a)
    ! Subroutine arguments
    real(c_double), dimension(:, :, :), target, intent(inout) :: a
    ! Local variables
    integer(4)  :: status,n1,n2,n3
    integer(8)  :: Nfft
    real(8)     :: bscale
    type(c_ptr) :: fft_ptr
    integer(4), dimension(3)              :: Lfft
    integer(4), dimension(4)              :: strides_in,strides_out
    type(dfti_descriptor), pointer        :: desc
    real(c_double), pointer, dimension(:) :: a_fft

    ! Get array dimensions
    n1 = size(a, dim=1) - 2  ! FFT padding
    n2 = size(a, dim=2)
    n3 = size(a, dim=3)

    ! Create pointer
    Nfft    = int(n1 + 2, kind=8) * int(n2, kind=8) * int(n3, kind=8)
    fft_ptr = c_loc(a)
    call c_f_pointer(fft_ptr, a_fft, [Nfft])

    ! Initialize
    Lfft        = (/ n1, n2, n3 /)
    strides_in  = (/ 0, 1, n1 / 2 + 1, (n1 / 2 + 1) * n2 /)
    strides_out = (/ 0, 1, n1 + 2, (n1 + 2) * n2 /)
    bscale      = 1D0 / (int(n1, kind=8) * int(n2, kind=8) * int(n3, kind=8))

    ! Initialize descriptor
    status = DftiCreateDescriptor(desc, DFTI_DOUBLE, DFTI_REAL, 3, Lfft)
    status = DftiSetValue(desc, DFTI_CONJUGATE_EVEN_STORAGE, &
         DFTI_COMPLEX_COMPLEX)
    status = DftiSetValue(desc, DFTI_INPUT_STRIDES, strides_in)
    status = DftiSetValue(desc, DFTI_OUTPUT_STRIDES, strides_out)
    status = DftiSetValue(desc, DFTI_BACKWARD_SCALE, bscale)
    status = DftiCommitDescriptor(desc)

    ! Compute transform
    status = DftiComputeBackward(desc, a_fft)

    ! Free descriptor
    status = DftiFreeDescriptor(desc)

    return
  end subroutine ifft_3d_d


  subroutine fft_2d_s(a)
    ! Subroutine arguments
    real(c_float), dimension(:, :), target, intent(inout) :: a
    ! Local variables
    integer(4)  :: status,n1,n2
    integer(8)  :: Nfft
    type(c_ptr) :: fft_ptr
    integer, dimension(2) :: Lfft
    integer, dimension(3) :: strides_in,strides_out
    type(dfti_descriptor), pointer :: desc
    real(c_float), pointer, dimension(:) :: a_fft

    ! Get array dimensions
    n1 = size(a, dim=1) - 2  ! FFT padding
    n2 = size(a, dim=2)

    ! Create pointer
    Nfft = int(n1 + 2, kind=8) * int(n2, kind=8)
    fft_ptr = c_loc(a)
    call c_f_pointer(fft_ptr, a_fft, [Nfft])

    ! Initialize
    Lfft        = (/ n1, n2 /)
    strides_in  = (/ 0, 1, n1 + 2 /)
    strides_out = (/ 0, 1, n1 / 2 + 1 /)

    ! Initialize descriptor
    status = DftiCreateDescriptor(desc, DFTI_SINGLE, DFTI_REAL, 2, Lfft)
    status = DftiSetValue(desc, DFTI_CONJUGATE_EVEN_STORAGE, &
         DFTI_COMPLEX_COMPLEX)
    status = DftiSetValue(desc, DFTI_INPUT_STRIDES, strides_in)
    status = DftiSetValue(desc, DFTI_OUTPUT_STRIDES, strides_out)
    status = DftiCommitDescriptor(desc)

    ! Compute transform
    status = DftiComputeForward(desc, a_fft)

    ! Free descriptor
    status = DftiFreeDescriptor(desc)

    return
  end subroutine fft_2d_s


  subroutine ifft_2d_s(a)
    ! Subroutine arguments
    real(c_float), dimension(:, :), target, intent(inout) :: a
    ! Local variables
    integer(4)  :: status,n1,n2
    integer(8)  :: Nfft
    real(8)     :: bscale
    type(c_ptr) :: fft_ptr
    integer, dimension(2) :: Lfft
    integer, dimension(3) :: strides_in,strides_out
    type(dfti_descriptor), pointer :: desc
    real(c_float), pointer, dimension(:) :: a_fft

    ! Get array dimensions
    n1 = size(a, dim=1) - 2  ! FFT padding
    n2 = size(a, dim=2)

    ! Create pointer
    Nfft = int(n1 + 2, kind=8) * int(n2, kind=8)
    fft_ptr = c_loc(a)
    call c_f_pointer(fft_ptr, a_fft, [Nfft])

    ! Initialize
    Lfft        = (/ n1, n2 /)
    strides_in  = (/ 0, 1, n1 / 2 + 1 /)
    strides_out = (/ 0, 1, n1 + 2 /)
    bscale = 1D0 / (int(n1, kind=8) * int(n2, kind=8))

    ! Initialize descriptor
    status = DftiCreateDescriptor(desc, DFTI_SINGLE, DFTI_REAL, 2, Lfft)
    status = DftiSetValue(desc, DFTI_CONJUGATE_EVEN_STORAGE, &
         DFTI_COMPLEX_COMPLEX)
    status = DftiSetValue(desc, DFTI_INPUT_STRIDES, strides_in)
    status = DftiSetValue(desc, DFTI_OUTPUT_STRIDES, strides_out)
    status = DftiSetValue(desc, DFTI_BACKWARD_SCALE, bscale)
    status = DftiCommitDescriptor(desc)

    ! Compute transform
    status = DftiComputeBackward(desc, a_fft)

    ! Free descriptor
    status = DftiFreeDescriptor(desc)

    return
  end subroutine ifft_2d_s


  subroutine fft_2d_d(a)
    ! Subroutine arguments
    real(c_double), dimension(:, :), target, intent(inout) :: a
    ! Local variables
    integer(4)  :: status,n1,n2
    integer(8)  :: Nfft
    type(c_ptr) :: fft_ptr
    integer, dimension(2)                 :: Lfft
    integer, dimension(3)                 :: strides_in,strides_out
    type(dfti_descriptor), pointer        :: desc
    real(c_double), pointer, dimension(:) :: a_fft

    ! Get array dimensions
    n1 = size(a, dim=1) - 2  ! FFT padding
    n2 = size(a, dim=2)

    ! Create pointer
    Nfft    = int(n1 + 2, kind=8) * int(n2, kind=8)
    fft_ptr = c_loc(a)
    call c_f_pointer(fft_ptr, a_fft, [Nfft])

    ! Initialize
    Lfft        = (/ n1, n2 /)
    strides_in  = (/ 0, 1, n1 + 2 /)
    strides_out = (/ 0, 1, n1 / 2 + 1 /)

    ! Initialize descriptor
    status = DftiCreateDescriptor(desc, DFTI_DOUBLE, DFTI_REAL, 2, Lfft)
    status = DftiSetValue(desc, DFTI_CONJUGATE_EVEN_STORAGE, &
         DFTI_COMPLEX_COMPLEX)
    status = DftiSetValue(desc, DFTI_INPUT_STRIDES, strides_in)
    status = DftiSetValue(desc, DFTI_OUTPUT_STRIDES, strides_out)
    status = DftiCommitDescriptor(desc)

    ! Compute transform
    status = DftiComputeForward(desc, a_fft)

    ! Free descriptor
    status = DftiFreeDescriptor(desc)

    return
  end subroutine fft_2d_d


  subroutine ifft_2d_d(a)
    ! Subroutine arguments
    real(c_double), dimension(:, :), target, intent(inout) :: a
    ! Local variables
    integer(4)  :: status,n1,n2
    integer(8)  :: Nfft
    real(8)     :: bscale
    type(c_ptr) :: fft_ptr
    integer, dimension(2)                 :: Lfft
    integer, dimension(3)                 :: strides_in,strides_out
    type(dfti_descriptor), pointer        :: desc
    real(c_double), pointer, dimension(:) :: a_fft

    ! Get array dimensions
    n1 = size(a, dim=1) - 2  ! FFT padding
    n2 = size(a, dim=2)

    ! Create pointer
    Nfft    = int(n1 + 2, kind=8) * int(n2, kind=8)
    fft_ptr = c_loc(a)
    call c_f_pointer(fft_ptr, a_fft, [Nfft])

    ! Initialize
    Lfft        = (/ n1, n2 /)
    strides_in  = (/ 0, 1, n1 / 2 + 1 /)
    strides_out = (/ 0, 1, n1 + 2 /)
    bscale      = 1D0 / (int(n1, kind=8) * int(n2, kind=8))

    ! Initialize descriptor
    status = DftiCreateDescriptor(desc, DFTI_DOUBLE, DFTI_REAL, 2, Lfft)
    status = DftiSetValue(desc, DFTI_CONJUGATE_EVEN_STORAGE, &
         DFTI_COMPLEX_COMPLEX)
    status = DftiSetValue(desc, DFTI_INPUT_STRIDES, strides_in)
    status = DftiSetValue(desc, DFTI_OUTPUT_STRIDES, strides_out)
    status = DftiSetValue(desc, DFTI_BACKWARD_SCALE, bscale)
    status = DftiCommitDescriptor(desc)

    ! Compute transform
    status = DftiComputeBackward(desc, a_fft)

    ! Free descriptor
    status = DftiFreeDescriptor(desc)

    return
  end subroutine ifft_2d_d


end module fft_tools
