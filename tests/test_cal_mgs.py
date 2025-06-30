import os
import sys
import time
import numpy as np
import h5py

# Ensure local package is importable
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Mulguisin.mgs_class import mulguisin_type


def run_notebook_example():
    """Replicate the cal_mgs.ipynb workflow."""
    data_path = os.path.join(os.path.dirname(__file__), "..", "test", "rand_gal_data_210218_sig10.hdf5")
    with h5py.File(data_path, "r") as f:
        x1 = np.array(f["x_rsd"])
        y1 = np.array(f["y_rsd"])
        z1 = np.array(f["z_rsd"])

    mgs_cls = mulguisin_type("voronoi")
    mgs = mgs_cls(10.0, x1, y1, z1)
    t0 = time.time()
    Nmgs, imgs, clg, clm, cng = mgs.get_mgs()
    elapsed = time.time() - t0
    return Nmgs, elapsed


def test_notebook_dataset():
    Nmgs, elapsed = run_notebook_example()
    print(f"Notebook dataset -> groups: {Nmgs}, time: {elapsed:.6f}s")
    assert Nmgs == 50

