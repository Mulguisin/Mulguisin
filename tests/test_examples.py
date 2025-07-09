import os
import sys
import time
import numpy as np

# Ensure local package is importable
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Mulguisin.mgs_class import mulguisin


def run_2d_example(n=1000):
    """Run a 2D clustering example with ``n`` points."""
    rng = np.random.default_rng(0)
    half = n // 2
    cluster1 = rng.normal(0.0, 0.05, size=(half, 2))
    cluster2 = rng.normal(1.0, 0.05, size=(n - half, 2))
    data = np.vstack((cluster1, cluster2))
    mgs = mulguisin(0.3, data[:, 0], data[:, 1], isort=np.arange(n))
    t0 = time.time()
    Nmgs, imgs, clg, clm, cng = mgs.get_mgs()
    elapsed = time.time() - t0
    return Nmgs, elapsed

def run_3d_example(n=1000):
    """Run a 3D clustering example with ``n`` points."""
    rng = np.random.default_rng(1)
    half = n // 2
    cluster1 = rng.normal(0.0, 0.05, size=(half, 3))
    cluster2 = rng.normal(1.0, 0.05, size=(n - half, 3))
    data = np.vstack((cluster1, cluster2))
    mgs = mulguisin(0.3, data[:, 0], data[:, 1], z1=data[:, 2], isort=np.arange(n))
    t0 = time.time()
    Nmgs, imgs, clg, clm, cng = mgs.get_mgs()
    elapsed = time.time() - t0
    return Nmgs, elapsed


def test_mgs_2d():
    Nmgs, elapsed = run_2d_example()
    print(f"2D example -> groups: {Nmgs}, time: {elapsed:.6f}s")
    assert Nmgs == 2


def test_mgs_3d():
    Nmgs, elapsed = run_3d_example()
    print(f"3D example -> groups: {Nmgs}, time: {elapsed:.6f}s")
    assert Nmgs == 2


def test_degree_stats():
    rng = np.random.default_rng(0)
    half = 50
    cluster1 = rng.normal(0.0, 0.05, size=(half, 2))
    cluster2 = rng.normal(1.0, 0.05, size=(half, 2))
    data = np.vstack((cluster1, cluster2))
    mgs = mulguisin(0.3, data[:, 0], data[:, 1], isort=np.arange(2 * half))
    mgs.get_mgs()
    avg, mx = mgs.Get_AllDegreeStats()
    assert np.allclose(avg, [1.96, 1.96])
    assert mx == [5, 6]
