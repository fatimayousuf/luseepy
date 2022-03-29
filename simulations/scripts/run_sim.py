#!/bin/env python
import numpy as np
from pathlib import Path
from pyuvsim import uvsim  # type: ignore
from simsetup import gen_uvbdict, gen_obsparams

def run(outpath):
    uvd = uvsim.run_uvsim(
        "sim_files/obsparam.yaml", return_uv=True, quiet=False
    )

    save_dict = {
        "data": uvd.data_array,
        "freq": uvd.freq_array,
        "lst": uvd.lst_array,
    }

    np.savez(outpath, **save_dict)


if __name__ == "__main__":
    from argparse import ArgumentParser

    #parser = ArgumentParser()
    #parser.add_argument("dir", type=str)  # path to uvbeam of antenna
    #args = parser.parse_args()
    #DIR = args.dir

    gen_uvbdict()
    gen_obsparams()
    run("results/test.npz")

        
