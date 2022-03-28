from astropy.coordinates import Galactic  # type: ignore
from astropy import units as u  # type: ignore
import healpy  # type: ignore
import numpy as np
from pygdsm import GlobalSkyModel2016  # type: ignore
import pyradiosky  # type: ignore
from typing import Optional
import warnings


def nside2npix(nside: int) -> int:
    return 12 * nside ** 2


def npix2nside(npix: int) -> int:
    nside = np.sqrt(npix / 12)
    notint = ~np.isclose(nside, int(nside))
    notpower = ~np.isclose(np.log2(nside), int(np.log2(nside)))
    if notint or notpower:
        raise ValueError(f"npix has an invalid value {npix}.")
    return int(nside)


def _degrade(sky_map: np.ndarray, nside_out: int) -> np.ndarray:
    return healpy.ud_grade(
        sky_map, nside_out=nside_out, order_in="RING", order_out="RING"
    )


class SkyMap:
    def __init__(
        self,
        frequency: float = 25,
        freq_unit: str = "MHz",
        nside: Optional[int] = 16,
        npix: Optional[int] = None,
        healpix_map: Optional[np.ndarray] = None,
        degrade: bool = True,
        base_name: str = "data/pygdsm",
    ):
        """
        Must at least specify one of nside, npix, or healpix_map. If they
        are incompatible, nside and npix will change.

        nside: int, nside of healpix map, must be power of 2
        npix: int, npix of healpix map, must be s.t. sqrt(npix/12)
        is power of 2
        healpix_map: np.ndarray, a sky map in RING healpix coords. Must have
        units of K.
        degrade: bool, whether to degrade the healpix map to the given nside
        """

        if healpix_map is not None:
            if degrade:
                if nside is None:
                    if npix is None:
                        warnings.warn(
                            "Cannot degrade map since nside and npix are not"
                            "specified.",
                            UserWarning,
                        )
                    else:
                        nside = npix2nside(npix)
                if nside is not None:
                    healpix_map = _degrade(healpix_map, nside)
            degrade = False
            npix = healpix_map.size
            nside = npix2nside(npix)
        elif nside is None:
            if npix is None:
                raise ValueError(
                    "Must specify at least one of nside, npix, and " "healpix_map."
                )
            else:
                nside = npix2nside(npix)
        elif npix is None:
            npix = nside2npix(nside)

        freq_conversion = {"Hz": 1e-6, "kHz": 1e-3, "MHz": 1.0, "GHz": 1e3}
        if freq_unit not in freq_conversion:
            raise ValueError(
                f"Invalid frequency unit {freq_unit}, must be in"
                f"{freq_conversion.keys()}"
            )
        frequency *= freq_conversion[freq_unit]
        freq_unit = "MHz"

        self.frequency = frequency
        self.freq_unit = freq_unit
        self.nside = nside
        self.npix = npix
        self.healpix_map = healpix_map
        self.degrade = degrade
        self.base_name = base_name

    def gen_gsm(self):
        gsm = GlobalSkyModel2016(freq_unit="MHz")
        gsm.generate(self.frequency)
        healpix_map = gsm.generated_map_data
        # degrading map
        if self.degrade:
            healpix_map = _degrade(healpix_map, self.nside)
        else:
            self.npix = len(healpix_map)
            self.nside = npix2nside(self.npix)

        self.healpix_map = healpix_map

    def make_pyradio_skymap(self):
        hpx_indices = np.arange(self.npix)
        stokes = u.Quantity(np.zeros((4, 1, self.npix)), unit=u.K)
        stokes[0, 0] = self.healpix_map * u.K  # U = Q=V=0 #XXX

        skymodel = pyradiosky.SkyModel(
            frame="galactic",
            stokes=stokes,
            component_type="healpix",
            spectral_type="spectral_index",
            hpx_inds=hpx_indices,
            hpx_order="ring",
            nside=self.nside,
            reference_frequency=self.frequency * np.ones(self.npix),
            spectral_index=-2.5 * np.ones(self.npix),
        )
        assert skymodel.check()
        skymodel.write_skyh5(self.base_name + f"_nside{self.nside}_ssi.h5")


if __name__ == "__main__":
    print("Making nside=16 fixed spectral index map...")
    sm = SkyMap(nside=16, degrade=True)
    if sm.healpix_map is None:
        sm.gen_gsm()
    sm.make_pyradio_skymap()
    print("Making full res fixed spectral index map...")
    sm = SkyMap(degrade=False)
    if sm.healpix_map is None:
        sm.gen_gsm()
    sm.make_pyradio_skymap()
