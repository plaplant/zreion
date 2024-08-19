module general_tools
  ! Default
  implicit none

  ! Parameters
  real(8), parameter :: pi = acos(-1D0)


contains


  pure function sinc(x)
    ! Function arguments
    real(8), intent(in) :: x
    real(8)             :: sinc

    if (abs(x) > 1D-6) then
       sinc = sin(x) / x
    else
       sinc = 1 - x**2 / 6
    endif

    return
  end function sinc


  pure function tophat(x)
    ! Function arguments
    real(8), intent(in) :: x
    real(8)             :: tophat

    if (abs(x) > 1d-6) then
       tophat = 3 * (sin(x) - cos(x) / x) / x**3
    else
       tophat = 1 - x**2 / 10
    endif

    return
  end function tophat


  pure function interpolate(x,f,a,b,flag)
    ! Function arguments
    integer(4),   intent(in) :: a,b
    real(8),      intent(in) :: x,f(:,:)
    character(3), intent(in) :: flag
    real(8)                  :: interpolate
    ! Local parameters
    real(8), parameter :: dxmin = epsilon(real(0.))
    ! Local variables
    integer(4) :: i,n,i1,i2
    real(8)    :: x1,x2,y1,y2,dx

    ! Get size of array
    n = size(f, dim=2)

    ! Function must be monotonic for binary search
    if (f(a,1) < f(a,n)) then
       i1 = 1
       i2 = n
       do while (i2-i1 > 1)
          i = (i1+i2)/2
          if (x > f(a,i)) then
             i1 = i
          else
             i2 = i
          endif
       enddo
    else
       i1 = 1
       i2 = n
       do while (i2-i1 > 1)
          i = (i1+i2)/2
          if (x < f(a,i)) then
             i1 = i
          else
             i2 = i
          endif
       enddo
    endif

    ! Interpolate in linear or logarithmic space
    if (flag == "lin") then
       x1 = f(a,i1)
       y1 = f(b,i1)
       x2 = f(a,i2)
       y2 = f(b,i2)
       dx = x2 - x1

       if (abs(dx) < dxmin) then
          interpolate = (y1+y2)/2
       else
          interpolate = y1+(y2-y1)*((x-x1)/(x2-x1))
       endif
    else if (flag == "log") then
       x1 = Dlog(f(a,i1))
       y1 = Dlog(f(b,i1))
       x2 = Dlog(f(a,i2))
       y2 = Dlog(f(b,i2))
       dx = x2 - x1

       if (abs(dx) < dxmin) then
          interpolate = exp((y1+y2)/2)
       else
          interpolate = exp(y1+(y2-y1)*((Dlog(x)-x1)/(x2-x1)))
       endif
    else
       interpolate = huge(0.)
    endif

    return
  end function interpolate


end module general_tools
