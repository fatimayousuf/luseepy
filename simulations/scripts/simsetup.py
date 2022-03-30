# from luseesky import __path__
from pathlib import Path
import yaml  # type: ignore


def gen_uvbdict(
    ant_path: str = None,
    telescope_coords: str = "(0., 180., 0.)",
    outpath: str = "sim_files/uvbeam.yaml"
):
    """
    Generate the uvbeam.yaml file.

    Parameters
        ant_path: str, path to the .uvbeam file of the antenna
        telescope_coords: str, moon coordinates (lon, lat, alt) of telescope
        outpath: str, path + filename of the uvbeam.yaml file that goes into
        pyuvsim
    """
    uvbdict = {
        "beam_paths": {
            0: "beams/data/lusee_004_N.fits",
            1: "beams/data/lusee_004_E.fits",
            2: "beams/data/lusee_004_S.fits",
            3: "beams/data/lusee_004_W.fits"
        },
        "freq_interp_kind": "linear",
        "telescope_location": f"{telescope_coords}",
        "world": "moon",
        "telescope_name": "LuSEE-Night",
        "x_orientation": "east",
    }
    with open(outpath, "w") as f:
        yaml.safe_dump(uvbdict, f)


def gen_obsparams(
    ant_model: str = None,
    nside: int = 64,
    outpath: str = "sim_files/obsparam.yaml",
):
    """
    Generate the obsparam dict.

    Parameters:
        ant_model: str, some defining characteristic of antenna. This defines
        the outfile name so picking something unique lessens the chances of
        overwriting another result!
        nside: int, nside of healpix map of sky
        outpath: str, path + filename of obsparam.yaml file that gets fed into
        pyuvsim.
    """
    obsparams = {
        "filing": {
            "outdir": f"outputs/filing",
            "outfile_name": f"whatisthis_uvsim",
            "output_format": "uvfits",
        },
        "freq": {
            "Nfreqs": 50,
            "start_freq": 1000000.0,
            "end_freq": 50000000.0,
            "channel_width": 1000000.0,
        },
        "sources": {
            "catalog": f"sky_models/data/const_1e4K_16.hdf5"
            },
        "telescope": {
            "array_layout": f"sim_files/layout.dat",
            "telescope_config_name": f"sim_files/uvbeam.yaml",
        },
        "time": {"Ntimes": 4, "start_time": 2459630.0, "duration_days": 28},
        "select": {"bls": "[(0, 0),(1,1),(2,2),(3,3),(1,2)]"},
    }
    with open(outpath, "w") as f:
        yaml.safe_dump(obsparams, f)

        
