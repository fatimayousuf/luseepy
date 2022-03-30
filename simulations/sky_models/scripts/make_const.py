from astropy.coordinates import Galactic  # type: ignore
from astropy import units as u  # type: ignore
import numpy as np
import healpy  # type: ignore
import pyradiosky  # type: ignore


def nside2npix(nside: int) -> int:
    return 12 * nside ** 2


def write_const_T_map(T0, fname, nside=16):
    npix = nside2npix(nside)
    hpx_indices = np.arange(npix)
    stokes = u.Quantity(np.zeros((4, 1, npix)), unit=u.K)
    stokes[0, 0] = T0 * u.K
    frequency = 10.0  # naturally in MHz
    skymodel = pyradiosky.SkyModel(
        stokes=stokes,
        component_type="healpix",
        spectral_type="spectral_index",
        hpx_inds=hpx_indices,
        hpx_order="ring",
        nside=nside,
        reference_frequency=frequency * np.ones(npix),
        spectral_index=-2.5 * np.ones(npix),
    )
    assert skymodel.check()
    skymodel.write_skyh5(fname)


if __name__ == "__main__":
    print("Making T=5*10^4k monopole only map...")
    write_const_T_map(5e4, "data/const_5e4K_16.hdf5", nside=16)
