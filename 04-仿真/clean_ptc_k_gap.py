"""Reproduce the clean two-layer PTC k-gap.

Change the parameters below, then run:

    python "04-仿真/clean_ptc_k_gap.py"
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from propagation_matrix import C_0, ETA_0, P_j


# Model parameters from Yang et al. Figure 1(c).
n_1 = 4.0  # Refractive index of temporal layer 1.
n_2 = 2.0  # Refractive index of temporal layer 2.
t_bar = 1e-15  # tau_j / n_j, in seconds
normalized_k_max = 5.2  # Upper bound of 2*k*c*t_bar/pi.
num_k = 2001  # Number of k samples.

output_path = Path(__file__).with_name("clean_ptc_k_gap.png")  # Saved figure.


def calculate_dispersion():
    """Calculate the Floquet dispersion from M(k) = P_2(k) @ P_1(k)."""
    # The chiral constraint fixes each layer duration.
    tau_1 = n_1 * t_bar
    tau_2 = n_2 * t_bar

    normalized_k_values = np.linspace(
        0.0, normalized_k_max, num_k
    )
    k_values = normalized_k_values * np.pi / (2.0 * C_0 * t_bar)
    cos_omega_t = np.empty(num_k)

    # Propagate through layers 1 and 2 in chronological order.
    for index, k in enumerate(k_values):
        P_1 = P_j(k, n_1, tau_1, C_0, ETA_0)
        P_2 = P_j(k, n_2, tau_2, C_0, ETA_0)
        period_matrix = P_2 @ P_1
        cos_omega_t[index] = np.real(np.trace(period_matrix) / 2.0)

    # Real Floquet frequency exists only where |Tr(M)/2| <= 1.
    allowed_mask = np.abs(cos_omega_t) <= 1.0
    k_gap_mask = ~allowed_mask

    normalized_frequency = np.full(num_k, np.nan)
    normalized_frequency[allowed_mask] = (
        np.arccos(np.clip(cos_omega_t[allowed_mask], -1.0, 1.0))
        / (2.0 * np.pi)
    )

    return (
        normalized_k_values,
        normalized_frequency,
        cos_omega_t,
        k_gap_mask,
    )


def plot_dispersion(
    normalized_k_values,
    normalized_frequency,
    cos_omega_t,
    k_gap_mask,
):
    """Plot the dispersion and the corresponding cos(Omega T)."""
    figure, (dispersion_axis, cosine_axis) = plt.subplots(
        1,
        2,
        figsize=(11, 4.5),
        sharex=True,
    )

    dispersion_axis.fill_between(
        normalized_k_values,
        -0.5,
        0.5,
        where=k_gap_mask,
        color="lightgray",
        label=r"$k$-gap",
    )
    dispersion_axis.plot(
        normalized_k_values,
        normalized_frequency,
        color="black",
        label="Allowed band",
    )
    dispersion_axis.plot(
        normalized_k_values,
        -normalized_frequency,
        color="black",
    )

    dispersion_axis.set_xlabel(r"$2kc\bar{t}/\pi$")
    dispersion_axis.set_ylabel(r"$\Omega T/(2\pi)$")
    dispersion_axis.set_xlim(0.0, normalized_k_max)
    dispersion_axis.set_ylim(-0.52, 0.52)
    dispersion_axis.set_title(
        rf"Clean PTC: $n_1={n_1:g}$, $n_2={n_2:g}$, "
        rf"$\bar{{t}}={t_bar/1e-15:g}$ fs"
    )
    dispersion_axis.legend()
    dispersion_axis.grid(alpha=0.3)

    cosine_axis.fill_between(
        normalized_k_values,
        np.min(cos_omega_t),
        np.max(cos_omega_t),
        where=k_gap_mask,
        color="lightgray",
    )
    cosine_axis.plot(normalized_k_values, cos_omega_t, color="black")
    cosine_axis.axhline(1.0, color="gray", linestyle="--")
    cosine_axis.axhline(-1.0, color="gray", linestyle="--")
    cosine_axis.set_xlabel(r"$2kc\bar{t}/\pi$")
    cosine_axis.set_ylabel(r"$\cos(\Omega T)$")
    cosine_axis.set_xlim(0.0, normalized_k_max)
    cosine_axis.set_title(r"$\cos(\Omega T)=\mathrm{Tr}(M)/2$")
    cosine_axis.grid(alpha=0.3)

    figure.tight_layout()
    figure.savefig(output_path, dpi=160)
    plt.close(figure)


def main():
    (
        normalized_k_values,
        normalized_frequency,
        cos_omega_t,
        k_gap_mask,
    ) = calculate_dispersion()
    plot_dispersion(
        normalized_k_values,
        normalized_frequency,
        cos_omega_t,
        k_gap_mask,
    )
    print(f"Figure: {output_path.resolve()}")


if __name__ == "__main__":
    main()
