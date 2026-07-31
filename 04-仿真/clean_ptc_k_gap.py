"""Clean two-layer PTC dispersion and k-gap.

Default parameters reproduce the structure of Yang et al. Figure 1(c).
Change the parameters in the section below, then run:

    python "04-仿真/clean_ptc_k_gap.py"
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from propagation_matrix import C0, ETA0, P_j


# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

N1 = 4.0
N2 = 2.0
T_BAR = 1e-15  # tau_j / n_j, in seconds
X_MAX = 5.2  # x = 2*k*c*T_BAR/pi
NUM_K = 2001

OUTPUT_PATH = Path(__file__).with_name("clean_ptc_k_gap.png")


def analytical_chi(theta1, theta2, n1, n2):
    """Analytical result for chi(k) = Tr[M(k)]/2."""
    contrast = 0.5 * (n1 / n2 + n2 / n1)
    return (
        np.cos(theta1) * np.cos(theta2)
        - contrast * np.sin(theta1) * np.sin(theta2)
    )


def calculate_spectrum(n1, n2, t_bar, x_max, num_k):
    """Calculate the analytical and matrix Floquet results."""
    tau1 = n1 * t_bar
    tau2 = n2 * t_bar
    period = tau1 + tau2

    x = np.linspace(0.0, x_max, num_k)
    k = x * np.pi / (2.0 * C0 * t_bar)

    theta1 = k * C0 * tau1 / n1
    theta2 = k * C0 * tau2 / n2
    chi_analytical = analytical_chi(theta1, theta2, n1, n2)

    chi_matrix = np.empty(num_k, dtype=complex)
    multipliers = np.empty((2, num_k), dtype=complex)
    determinants = np.empty(num_k, dtype=complex)

    for index, wavevector in enumerate(k):
        p1 = P_j(wavevector, n1, tau1, C0, ETA0)
        p2 = P_j(wavevector, n2, tau2, C0, ETA0)

        # Layer 1 acts first, so it is on the right.
        period_matrix = p2 @ p1

        chi_matrix[index] = np.trace(period_matrix) / 2.0
        multipliers[:, index] = np.linalg.eigvals(period_matrix)
        determinants[index] = np.linalg.det(period_matrix)

    # In an allowed band |chi| <= 1 and Omega is real.
    tolerance = 1e-12
    allowed = np.abs(chi_analytical) <= 1.0 + tolerance
    k_gap = ~allowed

    normalized_frequency = np.full(num_k, np.nan)
    normalized_frequency[allowed] = (
        np.arccos(
            np.clip(np.real(chi_matrix[allowed]), -1.0, 1.0)
        )
        / (2.0 * np.pi)
    )

    return {
        "x": x,
        "period": period,
        "chi_analytical": chi_analytical,
        "chi_matrix": chi_matrix,
        "multipliers": multipliers,
        "determinants": determinants,
        "allowed": allowed,
        "k_gap": k_gap,
        "frequency": normalized_frequency,
    }


def check_results(result):
    """Print a few direct numerical checks."""
    chi_from_eigenvalues = np.sum(result["multipliers"], axis=0) / 2.0

    errors = {
        "analytical trace vs matrix trace": np.max(
            np.abs(
                result["chi_analytical"] - result["chi_matrix"]
            )
        ),
        "analytical trace vs eigenvalues": np.max(
            np.abs(
                result["chi_analytical"] - chi_from_eigenvalues
            )
        ),
        "det(M) vs 1": np.max(
            np.abs(result["determinants"] - 1.0)
        ),
        "mu1*mu2 vs 1": np.max(
            np.abs(np.prod(result["multipliers"], axis=0) - 1.0)
        ),
    }

    print("Clean PTC checks")
    for name, error in errors.items():
        print(f"  {name}: {error:.3e}")

    if max(errors.values()) > 1e-10:
        raise RuntimeError("The clean PTC numerical checks failed.")


def plot_results(result, n1, n2, t_bar, output_path):
    """Draw the dispersion, trace comparison, and one simple error curve."""
    x = result["x"]
    frequency = result["frequency"]
    k_gap = result["k_gap"]
    trace_error = np.abs(
        result["chi_analytical"] - result["chi_matrix"]
    )

    figure, axes = plt.subplots(
        3,
        1,
        figsize=(8, 9),
        sharex=True,
    )

    # Top: Figure-1(c)-style dispersion.
    axes[0].fill_between(
        x,
        -0.5,
        0.5,
        where=k_gap,
        color="lightgray",
        label=r"$k$-gap (growth/decay)",
    )
    axes[0].plot(x, frequency, color="black", label="Allowed band")
    axes[0].plot(x, -frequency, color="black")
    axes[0].set_ylabel(r"$\Omega T/(2\pi)$")
    axes[0].set_ylim(-0.52, 0.52)
    axes[0].set_title(
        rf"Clean PTC: $n_1={n1:g}$, $n_2={n2:g}$, "
        rf"$\bar{{t}}={t_bar/1e-15:g}$ fs"
    )
    axes[0].legend()
    axes[0].grid(alpha=0.3)

    # Bottom: the two calculations should lie on top of each other.
    axes[1].fill_between(
        x,
        -1.4,
        1.1,
        where=k_gap,
        color="lightgray",
    )
    axes[1].plot(
        x,
        result["chi_analytical"],
        label="Analytical formula",
    )
    axes[1].plot(
        x,
        np.real(result["chi_matrix"]),
        "--",
        label=r"Matrix $M=P_2P_1$",
    )
    axes[1].axhline(1.0, color="black", linestyle=":", linewidth=1)
    axes[1].axhline(-1.0, color="black", linestyle=":", linewidth=1)
    axes[1].set_ylabel(r"$\chi(k)=\mathrm{Tr}(M)/2$")
    axes[1].set_ylim(-1.4, 1.1)
    axes[1].legend()
    axes[1].grid(alpha=0.3)

    # Bottom: one direct absolute-error curve.
    axes[2].plot(
        x[::5],
        trace_error[::5],
        ".",
        markersize=3,
    )
    axes[2].set_xlabel(r"$2kc\bar{t}/\pi$")
    axes[2].set_ylabel("Absolute error")
    axes[2].set_title(
        r"$|\chi_{\mathrm{analytical}}-\chi_{\mathrm{matrix}}|$"
    )
    axes[2].grid(alpha=0.3)
    axes[2].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    figure.tight_layout()
    figure.savefig(output_path, dpi=160)
    plt.close(figure)


def main():
    result = calculate_spectrum(N1, N2, T_BAR, X_MAX, NUM_K)
    check_results(result)
    plot_results(result, N1, N2, T_BAR, OUTPUT_PATH)
    print(f"  Figure: {OUTPUT_PATH.resolve()}")


if __name__ == "__main__":
    main()
