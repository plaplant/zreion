"""Module for cosmology-related calculations."""

import os
import warnings
import numpy as np
import astropy.cosmology
import camb

class CosmoParameters(object):
    """
    An object for storing cosmological and simulation parameters.

    The resulting object contains both cosmological parameters that define the
    simulation, as well as the computation parameters. For the cosmological
    parameters, we assume a default cosmology of Planck18, though the user can
    define/override with their own parameters if desired.

    Parameters
    ----------
    ndm : int
        The number of dark matter particles to use in the simulation per linear
        dimension. This means that the total number of particles in the
        simulation will be Ndm**3.
    ngrid : int
        The number of grid cells to use for the gridded output quantities (e.g.,
        21 cm brightness temperature). This should be less than or equal to
        `ndm` to avoid Poisson/shot-noise errors of grid sampling.
    box : float
        The size of one linear dimension of the simulation volume, in units of
        Mpc/h. This means the total volume is equal to box**3.
    zmean_zre : float
        The mean value of reionization, one of the model parameters of zreion.
        Note that in general, this may not be the midpoint of reionization, when
        the neutral fraction is 50% ionized, because the ionization history may
        be asymmetric. However, this value is generally close to the midpoint of
        reionization.
    alpha_zre : float
        The power law index in the zreion bias relation. Note that this controls
        the duration of reionization (along with kb), where larger values of
        alpha generally lead to longer reionziation scenarios.
    kb_zre : float
        The reference k-value in the zreion bias relation, in units of h/Mpc.
        Note that this controls the duration of reionization (along with alpha),
        where larger values of kb generally lead to longer reionization
        scenarios.
    ncpu : int, optional
        The number of CPUs/processors to use for the simulation, parallelized
        using OpenMP. If not specified, this will default to the number of CPUs
        available.
    cosmo : astropy.cosmology.Cosmology, optional
        An astropy cosmology object, from which the default cosmological
        parameters are taken. If `None`, then we default to Planck18.
    H0 : float, optional
        The Hubble constant, in units of km/s/Mpc.
    omegam : float, optional
        The density of matter (baryons + cold dark matter) at redshift 0,
        normalized by the critical density.
    omegab : float, optional
        The density of baryons at redshift 0, normalized by the critical
        density.
    omegac : float, optional
        The density of cold dark matter at redshift 0, normalized by the
        critical density.
    ns : float, optional
        The spectral index of scalar fluctuations.
    sigma8 : float, optional
        The amplitude of the density field, smoothed using a spherical tophat
        with a radius of 8 Mpc/h.
    m_nu : float, optional
        The sum of neutrino masses, in eV.
    b0_zre : float, optional
        The overall normalization of the zreion bias relation. The default is
        the inverse of the critical density of halo formation, \delta_c = 1.686.
    rsmooth_zre : float, optional
        The radius of a spherical tophat window used to smooth the zreion field,
        in units of Mpc/h. Default value is 1 Mpc/h.
    """
    def __init__(
        self,
        ndm,
        ngrid,
        box,
        zmean_zre,
        alpha_zre,
        kb_zre,
        ncpu=None,
        cosmo=None,
        H0=None,
        omegam=None,
        omegab=None,
        omegac=None,
        ns=None,
        sigma8=None,
        m_nu=None,
        b0_zre=1.0 / 1.686,
        rsmooth_zre=1.0,
    ):
        # set simulation parameters
        self.ndm = ndm
        if ngrid > ndm:
            warnings.warn(
                "`ngrid` should be less than or equal to `ndm` to avoid "
                "Poisson/shot-noise errors"
            )
        self.ngrid = ngrid
        self.box = box
        if ncpu is None:
            self.ncpu = os.cpu_count()
        else:
            self.ncpu = ncpu

        # set cosmology parameters
        if cosmo is None:
            # Default to Planck18 parameters
            cosmo = astropy.cosmology.Planck18
        elif not isinstance(cosmo, astropy.cosmology.Cosmology):
            raise ValueError(
                "`cosmo` parameter must be an astropy Cosmology class or subclass"
            )

        # Default to values from cosmo object, overriding when specified
        # Hubble constant
        if H0 is None:
            self.H0 = float(cosmo.H0.value)  # km/s/Mpc
        else:
            # check that H0 has the right units
            if H0 < 1 or H0 > 100:
                raise ValueError(
                    "The value of H0 should be in units of km/s/Mpc, and around 70."
                )
            self.H0 = H0

        # Relative density parameters
        # We allow for flexibility in specifying Omega_M vs. Omega_B vs. Omega_CDM,
        # but require that Omega_M = Omega_B + Omega_CDM
        if omegam is None:
            omegam = cosmo.Om0
        if omegab is None:
            omegab = cosmo.Ob0
        if omegac is None:
            omegac = cosmo.Odm0

        if not np.isclose(omegac + omegab, omegam):
            raise ValueError("omegam must equal omegab + omegac")
        else:
            # we're good
            self.omegam = omegam
            self.omegab = omegab
            self.omegac = omegac

        # we assume a flat universe and set Omega_Lambda accordingly
        self.omegal = 1.0 - self.omegam

        # we also assume that we have no evolution of the dark energy equation
        # of state
        self.wde = -1.0

        # we can infer Omega_R from the CMB temperature
        self.tcmb0 = cosmo.Tcmb0.value
        self.omegar = 4.48e-7 * (1 + 0.69) * self.tcmb0**4 / self.hubble0**2

        # power law index of scalar perturbations
        if ns is None:
            try:
                self.ns = cosmo.meta["n"]
            except (AttributeError, KeyError):
                raise ValueError(
                    "ns is not defined on cosmology object in expected way; "
                    "try defining it explicitly"
                )
        else:
            self.ns = ns

        # sigma8 -- degenerate with amplitude of power spectrum
        if sigma8 is None:
            try:
                self.sigma8 = cosmo.meta["sigma8"]
            except (AttributeError, KeyError):
                raise ValueError(
                    "sigma8 is not defined on cosmology object in expected way; "
                    "try defining it explicitly"
                )
        else:
            self.sigma8 = sigma8

        # m_nu -- sum of neutrino masses in eV
        if m_nu is None:
            try:
                self.m_nu = float(np.sum(cosmo.m_nu).to("eV").value)
            except (AttributeError, KeyError):
                raise ValueError(
                    "m_nu is not defined on cosmology object in expected way; "
                    "try defining it explicitly"
                )
        else:
            self.m_nu = m_nu

        # set placeholder attribute
        self.pk_lin = None

        # parameters defining zreion method
        self.b0_zre = b0_zre
        self.zmean_zre = zmean_zre
        self.alpha_zre = alpha_zre
        self.kb_zre = kb_zre
        self.rsmooth_zre = rsmooth_zre

    @property
    def hubble0(self):
        return self.H0 / 100

    def calculate_initial_ps(self):
        """
        Calcualte the initial power spectrum using CAMB.

        Use CAMB to calculate the power spectrum P(k) based on the cosmological
        parameters specified. The resulting power spectrum is saved on the
        object, which can be passed to Fortran for calculating the initial
        conditions. The resulting array has two columns: k (units: h/Mpc) and
        P(k) (units: (Mpc/h)**3). Note that the array is also in "Fortran order"
        (column-major), where the axes are reversed with respect to the
        traditional NumPy convention (row-major).

        After running, this function saves the output in the `pk_lin` attribute.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Raises
        ------
        ValueError
            This is raised if the renormalization of the power spectrum does not
            match the specified value of sigma8.
        """
        # set parameters for CAMB
        ombh2 = self.omegab * self.hubble0**2
        omch2 = self.omegac * self.hubble0**2
        As = 2e-9  # initial guess, to be renormalized later with sigma8
        kf = 2 * np.pi / self.box  # fundamental mode
        kmin = kf / 10
        kmax = max((kf * self.ngrid / 2) * 10, 2.0)

        # calculate initial pass, and renormalize
        pars = camb.set_params(
            H0=self.H0, ombh2=ombh2, omch2=omch2, As=As, ns=self.ns
        )
        pars.set_matter_power(redshifts=[0.0], kmax=kmax)
        results = camb.get_results(pars)
        s8_fid = results.get_sigma8_0()
        As_new = As * self.sigma8**2 / s8_fid**2
        pars.InitPower.set_params(As=As_new, ns=self.ns)

        # make sure result is consistent
        results = camb.get_results(pars)
        s8_new = results.get_sigma8_0()
        if not np.isclose(s8_new, self.sigma8):
            raise ValueError(
                f"power spectrum renormalization not correct; expected {sigma8}, "
                f"got {s8_new}"
            )

        # now get the power spectrum and save on object
        npoints = 2000
        kh, _, pk = results.get_matter_power_spectrum(minkh=kmin, maxkh=kmax, npoints=npoints)
        # we are going to pass this array into a Fortran function, so we want it
        # to be "Fortran contiguous" (i.e., column-major)
        pk_lin = np.asfortranarray([kh, pk[0, :]], dtype=np.float64)
        self.pk_lin = pk_lin

        return
