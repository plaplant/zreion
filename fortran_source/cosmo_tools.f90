module cosmo_tools
  ! Imports for C linking
  use iso_c_binding

  ! General tools
  use general_tools

  ! Default
  implicit none

  ! Local parameters
  integer(c_int), parameter :: N_scalefactor = 1000000
  real(c_double), parameter :: H0_cgs        = 3.24086D-18
  real(c_double), parameter :: Mpc2cm        = 3.08560D+24
  real(c_double), parameter :: rhoc_cgs      = 1.87890D-29

  ! Cosmology variables
  integer(c_int) :: Ndm,Ngrid,Ncpu,Ndomain
  real(c_double) :: omegam,omegal,omegab,omegar
  real(c_double) :: hubble0,sigma8,nsinit,wde,Tcmb0,box
  real(c_double) :: zmean_zre,b0_zre,kb_zre,alpha_zre,Rsmooth_zre
  real(c_double) :: x_unit,time_unit,Dfactor(2),vfactor(2)


contains


  pure function E_of_a(a)
    ! Assumes we have a flat universe (Omega_k = 0)
    ! Function arguments
    real(8), intent(in) :: a
    real(8)             :: E_of_a
    ! Local variables
    real(8) :: om,ol,or

    ! Assign values
    om  = omegam
    ol  = omegal
    or  = omegar

    E_of_a = sqrt(om / a**3 + or / a**4 + ol / a**(3 * (1 + wde)))
    return
  end function E_of_a


  pure function D_of_z(z)
    ! Function arguments
    real(8), intent(in) :: z
    real(8)             :: D_of_z(2)
    ! Local variables
    real(8) :: om,ol,or,aeq
    real(8) :: a,ah,aw,af,H,q,c0,c1
    real(8) :: da,D,Dh,dDda,dDdah

    ! Unpack cosmology
    om  = omegam
    ol  = omegal
    or  = omegar
    aeq = or / om

    ! Initialize
    af   = 1D0 / (1 + z)
    a    = max(aeq, 1D-6)
    D    = a + 2D0 / 3 * aeq
    dDda = 1

    ! Integrate differential equation until desired redshift
    do while (a < af)
       da = min(max(af - a, 1D-8), 1D-6)

       aw    = a**(3 * (1 + wde))
       H     = sqrt(om / a**3 + or / a**4 + ol / aw)
       q     = (om / a**3 + 2 * or / a**4 - 2 * ol / aw) / H**2 / 2
       c0    = -1.5D0 * om / H**2 / a**5
       c1    = (2 - q) / a
       Dh    = D + dDda * da / 2
       dDdah = dDda - (c1 * dDda + c0 * D) * da / 2
       ah    = a + da / 2

       aw    = ah**(3 * (1 + wde))
       H     = sqrt(om / ah ** 3 + or / ah**4 + ol / aw)
       q     = (om / ah**3 + 2 * or / ah**4 - 2 * ol / aw) / H**2 / 2
       c0    = -1.5 * om / H**2 / ah**5
       c1    = (2 - q) / ah
       D     = D + dDdah * da
       dDda  = dDda - (c1 * dDdah + c0 * Dh) * da
       a     = a + da
    enddo

    ! Save placeholder values
    D_of_z(1) = D
    D_of_z(2) = dDda * a / D

    ! Finish integrating to z = 0 to find overall normalization
    do while (a < 1)
       da    = min(max(1 - a, 1D-8), 1D-6)

       aw    = a**(3 * (1 + wde))
       H     = sqrt(om / a**3 + or / a**4 + ol / aw)
       q     = (om / a**3 + 2 * or / a**4 - 2 * ol / aw) / H**2 / 2
       c0    = -1.5D0 * om / H**2 / a**5
       c1    = (2 - q) / a
       Dh    = D + dDda * da / 2
       dDdah = dDda - (c1 * dDda + c0 * D) * da / 2
       ah    = a + da / 2

       aw    = ah**(3 * (1 + wde))
       H     = sqrt(om / ah**3 + or / ah**4 + ol / aw)
       q     = (om / ah**3 + 2 * or / ah**4 - 2 * ol / aw) / H**2 / 2
       c0    = -1.5D0 * om / H**2 / ah**5
       c1    = (2 - q) / a
       D     = D + dDdah * da
       dDda  = dDda - (c1 * dDdah + c0 * Dh) * da
       a     = a + da
    enddo

    D_of_z(1) = D_of_z(1) / D

    return
  end function D_of_z


  subroutine calc_cosmo_vars(z)
    ! Subroutine arguments
    real(8), intent(in) :: z
    ! Local variables
    real(8) :: om,ol,or,h0,Hubble
    real(8) :: a,hsq,oma,D(2),d1,d2,f1,f2

    ! Unpack cosmology
    om = omegam
    ol = omegal
    or = omegar
    h0 = hubble0

    ! Calculate auxiliary factors
    a      = 1D0 / (1 + z)
    hsq    = om / a**3 + ol / a**(3 * (1 + wde)) + or / a**4
    Hubble = H0_cgs * h0 * sqrt(hsq)
    oma    = om / a**3 / hsq**2

    ! Grid to physical values (cgs)
    x_unit    = (a * box / h0 / Ndm) * Mpc2cm
    time_unit = 2 / (3 * H0_cgs * h0 * sqrt(om)) * a**2

    ! Density and velocity factors
    D  = D_of_z(z)
    d1 = D(1)
    d2 = -3D0 / 7 * d1**2 * oma**(-1D0 / 143)
    f1 = D(2)
    f2 = 2 * oma**(6D0 / 11)
    Dfactor(1) = d1
    Dfactor(2) = d2
    vfactor(1) = d1 * f1 * Hubble
    vfactor(2) = d2 * f2 * Hubble

    return
  end subroutine calc_cosmo_vars


end module cosmo_tools
