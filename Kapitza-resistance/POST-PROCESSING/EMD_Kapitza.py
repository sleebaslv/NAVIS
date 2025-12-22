#!/usr/bin/env python3
"""
Kapitza resistance / interfacial thermal conductance analysis from LAMMPS outputs.


Expected input files (in --diranl):
  - tempdropdelta.<irun>.dat
  - fluxcnt.<irun>.dat

Output (in --results-dir):
  - kapitzas.<irun>.dat    columns: [Lk4, Rk4, avg]

Example (same positional args as your original):
  python kapitza_analysis.py 1 10 2000000 4 2000 5.0


"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Tuple

import numpy as np
from scipy import signal
from scipy.optimize import curve_fit


# ---------------------------- Helpers ---------------------------- #

def compute_area(diameter: float, length: float) -> float:
    """Compute interface area A = pi * d * L (matches original)."""
    return np.pi * diameter * length


def load_and_downsample(path: Path, dfreq: int, skiprows: int = 2) -> np.ndarray:
    """Load a 1D signal from text and downsample by taking every dfreq-th row."""
    data = np.loadtxt(path, skiprows=skiprows)
    return data[::dfreq]


def correlate_segmented(data1: np.ndarray, data2: np.ndarray, nrun: int) -> np.ndarray:
    """
    Split data into nrun equal segments, compute unbiased correlation per segment,
    average across segments. Returns length = segment_length.

    Unbiased normalization matches your original:
        corr_pos / counts, where counts = M, M-1, ..., 1 for segment length M.
    """
    if nrun <= 0:
        raise ValueError("nrun must be positive.")
    if data1.shape != data2.shape:
        raise ValueError("data1 and data2 must have the same shape.")
    if data1.ndim != 1:
        raise ValueError("data1 and data2 must be 1D arrays.")
    if len(data1) % nrun != 0:
        raise ValueError(
            f"len(data)={len(data1)} must be divisible by nrun={nrun} so splits are equal."
        )

    chunks1 = np.array_split(data1, nrun)
    chunks2 = np.array_split(data2, nrun)
    seg_len = len(chunks1[0])

    corr_all = np.zeros((nrun, seg_len), dtype=float)

    for j in range(nrun):
        corr_full = np.correlate(chunks1[j], chunks2[j], mode="full")
        n_full = corr_full.size

        # Positive lags (including lag 0)
        corr_pos = corr_full[n_full // 2 :]

        # Unbiased normalization by overlap count
        counts = np.arange(seg_len, 0, -1, dtype=float)  # M, M-1, ..., 1
        corr_unbiased = corr_pos / counts

        corr_all[j, :] = corr_unbiased

    return corr_all.mean(axis=0)


def moving_window_smooth(x: np.ndarray, window: int) -> np.ndarray:
    """
    Smooth x with a Hann window:
        convolve(x, hann, same) / sum(hann)
    """
    if window <= 1:
        return x.copy()
    win = signal.windows.hann(window)
    return signal.convolve(x, win, mode="same") / win.sum()


def trapz_integral(y: np.ndarray, dx: float) -> float:
    """Convenience wrapper for np.trapz."""
    return float(np.trapz(y, dx=dx))


# ---------------------------- Laplace + model fit ---------------------------- #

def laplace_transform(signal_t: np.ndarray, time: np.ndarray, s: np.ndarray) -> np.ndarray:
    r"""
    Numerical (finite-time) Laplace transform:

    \[
        F(s) = \int C(t) e^{-s t} dt
    \]

    This matches your original loop exactly (trapz of exp(-s*t) kernel).
    """
    out = np.zeros_like(s, dtype=float)
    for i in range(s.size):
        out[i] = np.trapz(signal_t * np.exp(-s[i] * time), time)
    return out


def kapitza_model(X: Tuple[np.ndarray, np.ndarray], B: float, L: float) -> np.ndarray:
    """
    Kapitza interfacial response model (your n=1 form):

        jt = (B * jj) / (s + L)

    where X = (jj, s).
    """
    jj, s = X
    return (B * jj) / (s + L)


# ---------------------------- Main analysis ---------------------------- #

def run_analysis(
    irun: str,
    dfreq: int,
    prodrun: int,
    nrun: int,
    window: int,
    diameter: float,
    diranl: Path,
    results_dir: Path,
    length_L: float = 51.578,
    temperature_K: float = 300.0,
    freq_max: float = 2000.0,
    freq_points: int = 200_000,
) -> np.ndarray:
    """
    Perform the full analysis. Returns array [[Lk4, Rk4, avg]] and writes output file.
    """
    # Derived values
    dt = 1.0 * dfreq
    area = compute_area(diameter, length_L)

    # After downsampling, expected number of points (matches original tsteps)
    tsteps = int(prodrun / dfreq)
    if tsteps <= 0:
        raise ValueError("Computed tsteps <= 0. Check prodrun and dfreq.")
    if tsteps % nrun != 0:
        raise ValueError(f"tsteps={tsteps} must be divisible by nrun={nrun}.")

    trim = int(tsteps / nrun)

    # Load data
    temp_path = diranl / f"tempdropdelta.{irun}.dat"
    flux_path = diranl / f"fluxcnt.{irun}.dat"

    temp_delta = load_and_downsample(temp_path, dfreq=dfreq, skiprows=2)
    flux_cnt = load_and_downsample(flux_path, dfreq=dfreq, skiprows=2)

    # Normalize flux by area (same as original)
    heat_flux = flux_cnt / area

    # Keep original behavior: L and R identical
    Ldelta = temp_delta
    Rdelta = temp_delta
    Lheatflux = heat_flux
    Rheatflux = heat_flux

    # Unit conversion constants (kept as-is)
    kb = 1.987204e-3  # real units
    kCal2Joule = 4184.0
    avaga = 6.02214179e23
    fs2s = 1.0e-15
    A2m = 1.0e-10
    K_real2SI = (A2m * A2m * avaga * fs2s) / kCal2Joule

    # Correlations
    LQ_tot = correlate_segmented(Ldelta, Lheatflux, nrun=nrun)[:trim]
    LQ_totemp = correlate_segmented(Ldelta, Ldelta, nrun=nrun)[:trim]
    RQ_tot = correlate_segmented(Rdelta, Rheatflux, nrun=nrun)[:trim]
    RQ_totemp = correlate_segmented(Rdelta, Rdelta, nrun=nrun)[:trim]

    # Smooth
    LQ_tot = moving_window_smooth(LQ_tot, window=window)
    RQ_tot = moving_window_smooth(RQ_tot, window=window)
    LQ_totemp = moving_window_smooth(LQ_totemp, window=window)
    RQ_totemp = moving_window_smooth(RQ_totemp, window=window)

    # Integrals (kept for parity with original; not all used later)
    _Lint_jj = trapz_integral(LQ_tot, dx=dt)
    _Rint_jj = trapz_integral(RQ_tot, dx=dt)
    _Lint_jt = trapz_integral(LQ_totemp, dx=dt)
    _Rint_jt = trapz_integral(RQ_totemp, dx=dt)

    # Normalized (kept)
    _LQ_tot_norm = LQ_tot / np.max(np.abs(LQ_tot))
    _RQ_tot_norm = RQ_tot / np.max(np.abs(RQ_tot))
    _LQ_totemp_norm = LQ_totemp / np.max(np.abs(LQ_totemp))
    _RQ_totemp_norm = RQ_totemp / np.max(np.abs(RQ_totemp))
    _Q_tot_norm = 0.5 * (_LQ_tot_norm + _RQ_tot_norm)
    _Q_totemp_norm = 0.5 * (_LQ_totemp_norm + _RQ_totemp_norm)

    # Green-Kubo (flux-flux)
    LQ_totnew = correlate_segmented(Lheatflux, Lheatflux, nrun=nrun)[:trim]
    RQ_totnew = correlate_segmented(Rheatflux, Rheatflux, nrun=nrun)[:trim]

    Lint_jjnew = trapz_integral(LQ_totnew, dx=dt)
    Rint_jjnew = trapz_integral(RQ_totnew, dx=dt)

    Lkg = (temperature_K * temperature_K * kb) / (Lint_jjnew * area)
    Rkg = (temperature_K * temperature_K * kb) / (Rint_jjnew * area)
    Lkg *= K_real2SI
    Rkg *= K_real2SI

    print(f"tsteps={tsteps} trim={trim} dt={dt}")
    print("Green-Kubo=", Rkg)
    print("Green-Kubo=", Lkg)

    # Btodd section: time grid must match trim length
    time = np.arange(0.0, dt * trim, dt)
    s = np.linspace(0.0, freq_max, num=freq_points)

    Lfs1 = laplace_transform(LQ_tot, time, s)
    Lfs2 = laplace_transform(LQ_totemp, time, s)
    Rfs1 = laplace_transform(RQ_tot, time, s)
    Rfs2 = laplace_transform(RQ_totemp, time, s)

    # Normalized (kept)
    _Lfs1_norm = Lfs1 / np.max(np.abs(Lfs1))
    _Lfs2_norm = Lfs2 / np.max(np.abs(Lfs2))
    _Rfs1_norm = Rfs1 / np.max(np.abs(Rfs1))
    _Rfs2_norm = Rfs2 / np.max(np.abs(Rfs2))

    # Fit kapitza_model
    guess = (1000.0, 1000.0)

    Lpopt, _ = curve_fit(kapitza_model, (Lfs1, s), Lfs2, p0=guess, maxfev=20000)
    Rpopt, _ = curve_fit(kapitza_model, (Rfs1, s), Rfs2, p0=guess, maxfev=20000)

    Lk4 = (Lpopt[0] / Lpopt[1]) * K_real2SI
    Rk4 = (Rpopt[0] / Rpopt[1]) * K_real2SI
    k4 = 0.5 * (Lk4 + Rk4)

    print("Kapitza_4=", Lk4)
    print("Kapitza_4=", Rk4)

    # Save output
    results_dir.mkdir(parents=True, exist_ok=True)
    out = np.array([[Lk4, Rk4, k4]], dtype=float)

    out_path = results_dir / f"kapitzas.{irun}.dat"
    np.savetxt(out_path, out)

    return out


# ---------------------------- CLI ---------------------------- #

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Compute Kapitza metrics from LAMMPS analysis files."
    )
    p.add_argument("irun", type=str, help="Run identifier used in input filenames.")
    p.add_argument("dfreq", type=int, help="Downsampling frequency (read every dfreq-th row).")
    p.add_argument("prodrun", type=int, help="Total production run length (original timestep units).")
    p.add_argument("nrun", type=int, help="Number of segments for correlation averaging.")
    p.add_argument("window", type=int, help="Hann window length for smoothing.")
    p.add_argument("dia", type=float, help="Diameter for area calculation (A = pi*d*L).")

    p.add_argument(
        "--diranl",
        type=Path,
        default=Path("../LAMMPS/"),
        help="Directory containing LAMMPS output files.",
    )
    p.add_argument(
        "--results-dir",
        type=Path,
        default=Path("./results/"),
        help="Directory to write kapitzas.<irun>.dat",
    )

    # Optional but useful knobs
    p.add_argument("--L", type=float, default=51.578, help="Length used for area calculation.")
    p.add_argument("--T", type=float, default=300.0, help="Temperature in K (Green-Kubo).")
    p.add_argument("--freq-max", type=float, default=2000.0, help="Max s-grid value for Laplace transform.")
    p.add_argument("--freq-points", type=int, default=200_000, help="Number of s-grid points.")
    return p


def main() -> None:
    args = build_parser().parse_args()

    run_analysis(
        irun=args.irun,
        dfreq=args.dfreq,
        prodrun=args.prodrun,
        nrun=args.nrun,
        window=args.window,
        diameter=args.dia,
        diranl=args.diranl,
        results_dir=args.results_dir,
        length_L=args.L,
        temperature_K=args.T,
        freq_max=args.freq_max,
        freq_points=args.freq_points,
    )


if __name__ == "__main__":
    main()
