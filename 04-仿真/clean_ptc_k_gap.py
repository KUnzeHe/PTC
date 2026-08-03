"""Reproduce the clean two-layer PTC k-gap.

Change the parameters below, then run:

    python "04-仿真/clean_ptc_k_gap.py"
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from propagation_matrix import C0, ETA0, P_j


# Parameters from Yang et al. Figure 1(c)
N1 = 4.0
N2 = 2.0
T_BAR = 1e-15  # tau_j / n_j, in seconds
X_MAX = 5.2  # x = 2*k*c*T_BAR/pi
NUM_K = 2001

OUTPUT_PATH = Path(__file__).with_name("clean_ptc_k_gap.png")


def calculate_dispersion():
    """Calculate the Floquet dispersion from M(k) = P2(k) @ P1(k)."""
    tau1 = N1 * T_BAR
    tau2 = N2 * T_BAR

    x = np.linspace(0.0, X_MAX, NUM_K)
    k = x * np.pi / (2.0 * C0 * T_BAR)
    cos_omega_T = np.empty(NUM_K)

    for index, wavevector in enumerate(k):
        p1 = P_j(wavevector, N1, tau1, C0, ETA0)
        p2 = P_j(wavevector, N2, tau2, C0, ETA0)
        period_matrix = p2 @ p1
        cos_omega_T[index] = np.real(np.trace(period_matrix) / 2.0)

    allowed = np.abs(cos_omega_T) <= 1.0
    k_gap = ~allowed

    frequency = np.full(NUM_K, np.nan)
    frequency[allowed] = (
        np.arccos(np.clip(cos_omega_T[allowed], -1.0, 1.0))
        / (2.0 * np.pi)
    )

    return x, frequency, cos_omega_T, k_gap


def plot_dispersion(x, frequency, cos_omega_T, k_gap):
    """Plot the dispersion and the corresponding cos(Omega T)."""
    figure, (dispersion_axis, cosine_axis) = plt.subplots(
        1,
        2,
        figsize=(11, 4.5),
        sharex=True,
    )

    dispersion_axis.fill_between(
        x,
        -0.5,
        0.5,
        where=k_gap,
        color="lightgray",
        label=r"$k$-gap",
    )
    dispersion_axis.plot(x, frequency, color="black", label="Allowed band")
    dispersion_axis.plot(x, -frequency, color="black")

    dispersion_axis.set_xlabel(r"$2kc\bar{t}/\pi$")
    dispersion_axis.set_ylabel(r"$\Omega T/(2\pi)$")
    dispersion_axis.set_xlim(0.0, X_MAX)
    dispersion_axis.set_ylim(-0.52, 0.52)
    dispersion_axis.set_title(
        rf"Clean PTC: $n_1={N1:g}$, $n_2={N2:g}$, "
        rf"$\bar{{t}}={T_BAR/1e-15:g}$ fs"
    )
    dispersion_axis.legend()
    dispersion_axis.grid(alpha=0.3)

    cosine_axis.fill_between(
        x,
        np.min(cos_omega_T),
        np.max(cos_omega_T),
        where=k_gap,
        color="lightgray",
    )
    cosine_axis.plot(x, cos_omega_T, color="black")
    cosine_axis.axhline(1.0, color="gray", linestyle="--")
    cosine_axis.axhline(-1.0, color="gray", linestyle="--")
    cosine_axis.set_xlabel(r"$2kc\bar{t}/\pi$")
    cosine_axis.set_ylabel(r"$\cos(\Omega T)$")
    cosine_axis.set_xlim(0.0, X_MAX)
    cosine_axis.set_title(r"$\cos(\Omega T)=\mathrm{Tr}(M)/2$")
    cosine_axis.grid(alpha=0.3)

    figure.tight_layout()
    figure.savefig(OUTPUT_PATH, dpi=160)
    plt.close(figure)


def main():
    x, frequency, cos_omega_T, k_gap = calculate_dispersion()
    plot_dispersion(x, frequency, cos_omega_T, k_gap)
    print(f"Figure: {OUTPUT_PATH.resolve()}")


if __name__ == "__main__":
    main()
