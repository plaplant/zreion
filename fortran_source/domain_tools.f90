module domain_tools
  ! OpenMP
  use omp_lib

  ! Cosmology tools
  use cosmo_tools

  ! Default
  implicit none

  ! Module arrays
  integer(4), allocatable, target, dimension(:,:)   :: domain
  integer(4), allocatable, target, dimension(:,:,:) :: indx_dom


contains


  subroutine init_domain
    ! Local variables
    integer(4) :: a,b,c
    integer(4) :: i,j,k
    integer(4) :: l,m,n
    integer(4) :: i1,j1,k1
    integer(4) :: n1,n2,n3,nx,ny,nz

    ! Calculate domain subdivision
    n = nint((Ncpu)**(1.0D0 / 3))

    do i=n,Ngrid
       if (mod(Ngrid, 2 * i) == 0) then
          n1 = 2 * i
          n2 = n1
          n3 = n1
          nx = Ngrid / n1
          ny = nx
          nz = nx
          Ndomain = n1 * n2 * n3
          exit
       endif
    enddo

    ! Allocate domain arrays
    if (allocated(domain)) then
       deallocate(domain)
       deallocate(indx_dom)
    endif
    allocate(domain(27, Ndomain))
    allocate(indx_dom(2, 3, Ndomain))

    ! Calculate domain information
    do n=1,Ndomain
       i = 1 + mod((n - 1), n1)
       j = 1 + mod((n - 1) / n1, n2)
       k = 1 +    ((n - 1) / n1 / n2)

       ! Calculate indices
       indx_dom(1, 1, n) = 1 + (i - 1) * nx
       indx_dom(2, 1, n) = i * nx
       indx_dom(1, 2, n) = 1 + (j - 1) * ny
       indx_dom(2, 2, n) = j * ny
       indx_dom(1, 3, n) = 1 + (k - 1) * nz
       indx_dom(2, 3, n) = k * nz

       ! Init domain
       m            = 1
       domain(1, n) = 0

       ! Calculate neighbors
       do c=-1,1
          k1 = 1 + mod(k + c - 1 + n3, n3)
          do b=-1,1
             j1 = 1 + mod(j + b - 1 + n2, n2)
             do a=-1,1
                i1 = 1 + mod(i + a - 1 + n1, n1)
                l  = i1 + (j1 - 1) * n1 + (k1 - 1) * n1 * n2

                if (l /= n) then
                   m            = m + 1
                   domain(m, n) = l
                endif
             enddo
          enddo
       enddo
    enddo

    return
  end subroutine init_domain


  subroutine set_domain(i, Ndomain, domain)
    ! Subroutine arguments
    integer(4), intent(inout) :: i
    integer(4), intent(in)    :: Ndomain
    integer(4), intent(inout) :: domain(:,:)
    ! Local variables
    integer(4) :: l,m,n
    logical    :: loop,eval

    !$omp critical (access_domain)
    ! Repeat until eligible domain is found
    loop = .true.

    do while (loop)
       do n=1,Ndomain
          if (domain(1, n) == 0) then
             eval = .true.
             do m=2,27
                l = domain(m, n)
                if (domain(1, l) == 1) then
                   eval = .false.
                   exit
                endif
             enddo

             if (eval) then
                i            = n
                domain(1, i) = 1
                loop         = .false.
                exit
             endif
          endif
       enddo
    enddo
    !$omp end critical (access_domain)

    return
  end subroutine set_domain


  subroutine end_domain(i, Ndomain, domain)
    ! Subroutine arguments
    integer(4), intent(in)    :: i,Ndomain
    integer(4), intent(inout) :: domain(:,:)

    ! Domain evaluated
    domain(1,i) = 2

    ! Reset domains if all are evaluated
    if (sum(domain(1, :)) == 2 * Ndomain) domain(1, :) = 0

    return
  end subroutine end_domain

end module domain_tools
