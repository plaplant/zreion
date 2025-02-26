"""Main module for interacting with simulation methods."""

import h5py
import numpy as np

from .cosmology import CosmoParameters
from . import fortran_tools


def set_cosmo_parameters(cosmo_params):
    """
    Function to set cosmological parameters in the Fortran module.

    Before performing any of the cosmology-related calculations in the Fortran
    module, we need to update the parameters accordingly. They are accessible as
    attributes in the `cosmo_tools` module, so they are easy to access and
    update.

    Parameters
    ----------
    cosmo_params : CosmoParameters object
        The cosmological and simulation parameters of the simulation.

    Returns
    -------
    None
    """
    # Update Fortran module-level variables via related attributes
    # Simulation parameters
    fortran_tools.cosmo_tools.ndm = cosmo_params.ndm
    fortran_tools.cosmo_tools.ngrid = cosmo_params.ngrid
    fortran_tools.cosmo_tools.ncpu = cosmo_params.ncpu

    # Cosmology parameters
    fortran_tools.cosmo_tools.omegam = cosmo_params.omegam
    fortran_tools.cosmo_tools.omegal = cosmo_params.omegal
    fortran_tools.cosmo_tools.omegab = cosmo_params.omegab
    fortran_tools.cosmo_tools.omegar = cosmo_params.omegar
    fortran_tools.cosmo_tools.hubble0 = cosmo_params.hubble0
    fortran_tools.cosmo_tools.sigma8 = cosmo_params.sigma8
    fortran_tools.cosmo_tools.nsinit = cosmo_params.ns
    fortran_tools.cosmo_tools.wde = cosmo_params.wde
    fortran_tools.cosmo_tools.tcmb0 = cosmo_params.tcmb0

    # Simulation parameters
    fortran_tools.cosmo_tools.box = cosmo_params.box
    fortran_tools.cosmo_tools.zmean_zre = cosmo_params.zmean_zre
    fortran_tools.cosmo_tools.b0_zre = cosmo_params.b0_zre
    fortran_tools.cosmo_tools.kb_zre = cosmo_params.kb_zre
    fortran_tools.cosmo_tools.alpha_zre = cosmo_params.alpha_zre
    fortran_tools.cosmo_tools.rsmooth_zre = cosmo_params.rsmooth_zre

    return


def generate_initial_conditions(cosmo_params, seed=None):
    """
    Generate initial conditions given simulation parameters.

    This function will use the specified cosmological parameters to generate the
    initial conditions suitable for a cosmological simulation. It does this by
    using CAMB to generate a (matter) transfer function given the cosmological
    parameters. It then uses this to generate the initial random phases that
    match the specified power spectrum. The resulting initial conditions are
    a series of packed Fourier coefficients using the "conjugate even storage"
    packing scheme. They are also in "Fortran" ordering, where the first index
    of the 3d array is the fastest-running index. Most of these issues should be
    handled automatically by NumPy and the underlying Fortran library.

    Parameters
    ----------
    cosmo_params : CosmoParameters object
        The cosmological and simulation parameters of the simulation.
    seed : int, optional
        The initial seed to use for generating random numbers. Allows for
        reproducible initial conditions when using the same seed. If not
        specified, a random seed is chosen between the numbers 1 (inclusive)
        and 2**31 (exclusive).

    Returns
    -------
    ics : ndarray of float
        The Fourier transformed initial conditions. This array should have shape
        (ndm + 2, ndm, ndm), where `ndm` is defined in CosmoParameters. The
        extra factor of 2 is from the conjugate-even storage used for saving the
        result of the Fourier transform.
    """
    # calculate initial power spectrum if necessary
    if cosmo_params.pk_lin is None:
        cosmo_params.calculate_initial_ps()

    # make an array for initial conditions
    ndm = cosmo_params.ndm
    ic_field = np.empty((ndm + 2, ndm, ndm), dtype=np.float64, order="F")

    # make a seed if the user didn't give us one
    if seed is None:
        rng = np.random.default_rng()
        seed = np.int32(rng.integers(1, 2**31))

    fortran_tools.ic_tools.generate_ics(ic_field, cosmo_params.pk_lin, np.int32(seed))

    return ic_field


def calc_density_and_velocity(redshift, cosmo_params, ic_field, dep_scheme=1):
    r"""
    Calculate the density and velocity fields at a given redshift.

    This function takes a particular redshift and a series of initial conditions
    and produces the corresponding density and velocity fields.

    Parameters
    ----------
    redshift : float
        The target redshift.
    cosmo_params : CosmoParameters object
        The cosmological and simulation parameters of the simulation.
    ic_field : ndarray of float
        The (Fourier transform) of the initial conditions of the simulation.
        Should have shape (ndm + 2, ndm, ndm), where `ndm` is defined in
        CosmoParameters.
    dep_scheme : int, optional
        The particle deposition scheme used for generating the density and
        velocity fields. Should be one of: 0 (for nearest grid point, NGP), 1
        (for cloud-in-cell, CIC), or 2 (for triangular-shaped clouds, TSC).

    Returns
    -------
    density : ndarray of float
        The density field, with shape (ngrid, ngrid, ngrid). This is the total
        matter density normalized by average cosmic density, :math:`1 + \delta`,
        where :math:`\delta` is the fractional overdensity. This field should
        have a mean of 1 and a minimum value of 0.
    velocity : ndarray of float
        The velocity field, which shape (3, ngrid, ngrid, ngrid). This is the
        proper velocity in units of km/s.
    """
    # check inputs
    redshift = float(redshift)
    if redshift < 0:
        raise ValueError("redshift must be positive")
    ic_ndm = ic_field.shape[0] - 2
    ndm = cosmo_params.ndm
    if ic_ndm != ndm:
        raise ValueError(
            f"ic_field is not the correct shape; expected {ndm}, got {ic_ndm}"
        )
    dep_scheme = int(dep_scheme)
    if dep_scheme not in (0, 1, 2):
        raise ValueError("dep_scheme must be one of: 0, 1, 2")

    # make arrays for density and velocity
    ngrid = cosmo_params.ngrid
    density = np.empty((ngrid, ngrid, ngrid), dtype=np.float64, order="F")
    velocity = np.empty((3, ngrid, ngrid, ngrid), dtype=np.float64, order="F")

    # call Fortran subroutine to compute fields
    fortran_tools.field_tools.calc_density_velocity_fields(
        redshift, ic_field, density, velocity, dep_scheme
    )

    return density, velocity


def calc_zreion(density, cosmo_params):
    """
    Calcualte the redshift of reionization field for a given simulation.

    This function uses the density field, calculated at the midpoint of
    reionization, and computes the resulting "redshift of reionization" field,
    which corresponds to the redshift at which that point in the volume was
    reionized. This field is needed to compute the ionization field or 21 cm
    field.

    Parameters
    ----------
    density : ndarray of float
        The density field of the midpoint of reionization.
    cosmo_params : CosmoParameters object
        The cosmological and simulation parameters of the simulation.

    Returns
    -------
    zreion : ndarray of float
        The redshift of reionization field.
    """
    # make arrays for zreion field
    # we need extra elements for FFT padding
    ngrid = cosmo_params.ngrid
    zreion = np.empty((ngrid + 2, ngrid, ngrid), dtype=np.float64, order="F")

    fortran_tools.zreion_tools.calc_zreion(density, zreion)

    return zreion


def calc_ion_t21(redshift, density, zreion, cosmo_params):
    """
    Calculate the ionization field and 21 cm brightness field.

    This function will calculate the ionization and 21 cm brightness fields at a
    given redshift.

    Parameters
    ----------
    redshift : float
        The target redshift.
    density : ndarray of float
        The density field of the midpoint of reionization.
    zreion : ndarray of float
        The redshift of reionization field.
    cosmo_params : CosmoParameters object
        The cosmological and simulation parameters of the simulation.

    Returns
    -------
    ion : ndarray of float
        The ionization field of the simulation at the specified redshift. Values
        of this field are 0 (for neutral) or 1 (for ionized).
    t21 : ndarray of float
        The 21 cm brightness temperature of the simulation at the specified
        redshift. Values are in millikelvin (mK).
    """
    # check input
    redshift = float(redshift)
    if redshift < 0:
        raise ValueError("redshift must be positive")
    ngrid = cosmo_params.ngrid
    density_ngrid = density.shape[0]
    if density_ngrid != ngrid:
        raise ValueError(
            f"density is not the correct shape; expected {ngrid}, got {density_ngrid}"
        )
    zreion_ngrid = zreion.shape[0] - 2
    if zreion_ngrid != ngrid:
        raise ValueError(
            f"zreion is not the correct shape; expected {ngrid}, got {zreion_ngrid}"
        )

    # make arrays for ionization and T21 fields
    ion = np.empty((ngrid, ngrid, ngrid), dtype=np.float64, order="F")
    t21 = np.empty((ngrid, ngrid, ngrid), dtype=np.float64, order="F")

    fortran_tools.field_tools.calc_ion_t21_fields(
        redshift, density, zreion, ion, t21
    )

    return ion, t21


class Simulation(object):
    r"""
    Convenience object for defining, running, and saving a simulation.

    This object has top-level methods for generating 21 cm brightness
    temperature fields at particular redshift values, as well as computing
    intermediate quantities that might be useful.

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
        the inverse of the critical density of halo formation,
        :math:`\delta_c = 1.686`.
    rsmooth_zre : float, optional
        The radius of a spherical tophat window used to smooth the zreion field,
        in units of Mpc/h. Default value is 1 Mpc/h.
    seed : int, optional
        The initial random seed to use for the simulation. If not specified,
        will default to a random number between 1 (inclusive) and 2**31
        (exclusive).
    dep_scheme : int, optional
        The particle deposition scheme used for generating the density and
        velocity fields. Should be one of: 0 (for nearest grid point, NGP), 1
        (for cloud-in-cell, CIC), or 2 (for triangular-shaped clouds, TSC).
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
        b0_zre=None,
        rsmooth_zre=None,
        seed=None,
        dep_scheme=None,
    ):
        # define a cosmology object and save it on the simulation
        cosmo_params = CosmoParameters(
            ndm=ndm,
            ngrid=ngrid,
            box=box,
            zmean_zre=zmean_zre,
            alpha_zre=alpha_zre,
            kb_zre=kb_zre,
            ncpu=ncpu,
            cosmo=cosmo,
            H0=H0,
            omegam=omegam,
            omegab=omegab,
            omegac=omegac,
            ns=ns,
            sigma8=sigma8,
            m_nu=m_nu,
            b0_zre=b0_zre,
            rsmooth_zre=rsmooth_zre,
        )
        self.cosmo_params = cosmo_params

        # define initial random seed
        if seed is None:
            rng = np.random.default_rng()
            seed = np.int32(rng.integers(1, 2**31))
        self.seed = seed

        # define self-consistent particle deposition scheme
        if dep_scheme is None:
            dep_scheme = 1
        dep_scheme = int(dep_scheme)
        if dep_scheme not in (0, 1, 2):
            raise ValueError("dep_scheme must be one of: 0, 1, 2")
        self.dep_scheme = dep_scheme

        # define placeholder attributes for later fields
        self.ic_field = None
        self.zreion = None
        self.redshift = None
        self.density = None
        self.velocity = None
        self.ion = None
        self.t21 = None

        return

    def generate_ics(self):
        """
        Generate initial conditions and save on the object.

        This function will generate initial conditions and save them on the
        object. Note that these are the Fourier-transformed initial
        fluctuations.

        Parameters
        ----------
        seed : int, optional
            The random seed to use for the initial conditions.

        Returns
        -------
        None
        """
        # set the cosmological parameters
        set_cosmo_parameters(self.cosmo_params)

        if self.ic_field is None:
            # generate the initial conditions
            self.ic_field = generate_initial_conditions(
                self.cosmo_params, seed=self.seed
            )

        return

    def calc_density_velocity(self, redshift):
        """
        Calculate the density and velocity fields at the given redshift.

        The resulting fields will be saved on the object.

        Parameters
        ----------
        redshift : float
            The redshift value to calculate the density and velocity fields at.

        Returns
        -------
        None
        """
        # validate input
        redshift = float(redshift)
        if redshift < 0:
            raise ValueError("redshift must be greater than or equal to 0")

        # generate initial conditions if needed
        if self.ic_field is None:
            self.generate_ics()

        self.redshift = redshift
        self.density, self.velocity = calc_density_and_velocity(
            redshift=self.redshift,
            cosmo_params=self.cosmo_params,
            ic_field=self.ic_field,
            dep_scheme=self.dep_scheme,
        )

        return

    def calc_zreion(self):
        """
        Calculate the redshift of reionization field zreion.

        If necessary, this function will first calculate the initial conditions
        field and the density field at the redshift specified by the zmean_zre
        parameter. These will be saved on the object for future use. This
        function will also save the `zreion` field on the object.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        # set the cosmological parameters
        set_cosmo_parameters(self.cosmo_params)

        # generate initial conditions if needed
        if self.ic_field is None:
            self.generate_ics()

        # compute density field at midpoint if needed
        if (self.redshift is None or
            not np.isclose(self.redshift, self.zmean_zre)):
            self.calc_density_velocity(self.zmean_zre)

        self.zreion = calc_zreion(self.density, self.cosmo_params)

        return

    def calc_ion_t21(self, redshift):
        """
        Calculate the ionization and cosmological 21 cm fields at some redshift.

        Parameters
        ----------
        redshift : float
            The redshift at which to calculate the ionization and 21 cm fields.

        Returns
        -------
        None
        """
        # validate input
        redshift = float(redshift)
        if redshift < 0:
            raise ValueError("redshift must be greater than or equal to 0")

        # set the cosmological parameters
        set_cosmo_parameters(self.cosmo_params)

        # generate initial conditions if needed
        if self.ic_field is None:
            self.generate_ics()

        # compute zreion field if needed
        if self.zreion is None:
            self.calc_zreion()

        # compute the requested density field
        if (self.redshift is None or
            not np.isclose(self.redshift, zmean_zre)):
            self.calc_density_velocity(zmean_zre)

        # compute the resulting fields and save on object
        self.ion, self.t21 = calc_ion_t21(
            redshift, self.density, self.zreion, self.cosmo_params
        )

        return
