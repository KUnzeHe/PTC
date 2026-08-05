"""Reproduce the clean chiral-PTC to SSH equivalence.

The script performs three independent checks for the unit-cell convention

    A_n --v-- B_n --w-- A_{n+1}.

It compares the TMM Floquet dispersion with the SSH bands, evaluates the
winding of z(q) = v + w exp(iq), and diagonalizes a finite open SSH chain.

Run from the project root with

    python "04-仿真/clean_ptc_ssh_equivalence.py"
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from propagation_matrix import C0, ETA0, P_j


# Yang et al. Figure 1(c): n1 > n2 is topological in the fixed convention.
REFERENCE_N1 = 4.0
REFERENCE_N2 = 2.0
T_BAR = 1e-15

NUM_K = 2001
NUM_Q = 1601
NUM_CELLS = 30

OUTPUT_FIGURE = Path(__file__).with_name(
    "clean_ptc_ssh_equivalence.png"
)


def _require_positive_finite(name, value):
    """Validate one positive scalar model parameter."""
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be positive and finite.")


def ssh_couplings(n1, n2):
    """Return fixed-convention intracell v and intercell w couplings."""
    _require_positive_finite("n1", n1)
    _require_positive_finite("n2", n2)
    denominator = n1 + n2
    return n2 / denominator, n1 / denominator


def ssh_bloch_hamiltonian(q, v, w):
    """Return the 2x2 clean SSH Bloch Hamiltonian at momentum q."""
    off_diagonal = v + w * np.exp(-1j * q)
    return np.array(
        [[0.0, off_diagonal], [np.conjugate(off_diagonal), 0.0]],
        dtype=complex,
    )


def ssh_bands(q, v, w):
    """Return the analytical lower and upper SSH bands."""
    band_magnitude = np.sqrt(v**2 + w**2 + 2.0 * v * w * np.cos(q))
    return -band_magnitude, band_magnitude


def calculate_tmm_mapping(n1, n2, t_bar, num_k=NUM_K):
    """Map TMM data from k to the common SSH coordinates (q, lambda).

    One theta interval [0, pi] is enough to cover both signs of lambda.
    In a k-gap, q is complex; those points are kept for the algebraic check
    but are excluded from the real-q overlay.
    """
    _require_positive_finite("t_bar", t_bar)
    if num_k < 3:
        raise ValueError("num_k must be at least 3.")

    tau1 = n1 * t_bar
    tau2 = n2 * t_bar
    theta = np.linspace(0.0, np.pi, num_k)
    wavevectors = theta / (C0 * t_bar)
    trace_half = np.empty_like(theta)

    for index, wavevector in enumerate(wavevectors):
        period_matrix = P_j(wavevector, n2, tau2, C0, ETA0) @ P_j(
            wavevector, n1, tau1, C0, ETA0
        )
        trace_half[index] = float(np.real(np.trace(period_matrix) / 2.0))

    lambda_from_k = np.cos(theta)
    tolerance = 50.0 * np.finfo(float).eps
    allowed = np.abs(trace_half) <= 1.0 + tolerance
    real_q = np.full_like(theta, np.nan)
    real_q[allowed] = np.arccos(np.clip(trace_half[allowed], -1.0, 1.0))

    v, w = ssh_couplings(n1, n2)
    trace_half_from_ssh = (
        lambda_from_k**2 - v**2 - w**2
    ) / (2.0 * v * w)

    return {
        "theta": theta,
        "wavevectors": wavevectors,
        "lambda": lambda_from_k,
        "trace_half": trace_half,
        "trace_half_from_ssh": trace_half_from_ssh,
        "allowed": allowed,
        "q": real_q,
    }


def winding_number(v, w, num_q=NUM_Q):
    """Compute the winding of z(q)=v+w*exp(iq) by phase increments."""
    if num_q < 5:
        raise ValueError("num_q must be at least 5.")
    if np.isclose(v, w, rtol=0.0, atol=100.0 * np.finfo(float).eps):
        raise ValueError("The winding number is undefined at the gap closing v=w.")

    q = np.linspace(0.0, 2.0 * np.pi, num_q)
    path = v + w * np.exp(1j * q)
    phase_increments = np.angle(path[1:] / path[:-1])
    winding = float(np.sum(phase_increments) / (2.0 * np.pi))
    return q, path, winding


def build_open_ssh_hamiltonian(num_cells, v, w):
    """Build a finite SSH chain with no B_N to A_1 closing bond."""
    if not isinstance(num_cells, int) or num_cells < 2:
        raise ValueError("num_cells must be an integer of at least 2.")

    hamiltonian = np.zeros((2 * num_cells, 2 * num_cells), dtype=float)
    for cell in range(num_cells):
        a_site = 2 * cell
        b_site = a_site + 1
        hamiltonian[a_site, b_site] = v
        hamiltonian[b_site, a_site] = v

        if cell < num_cells - 1:
            next_a_site = a_site + 2
            hamiltonian[b_site, next_a_site] = w
            hamiltonian[next_a_site, b_site] = w

    return hamiltonian


def open_chain_data(num_cells, v, w):
    """Return spectrum and the combined density of the closest-to-zero pair."""
    hamiltonian = build_open_ssh_hamiltonian(num_cells, v, w)
    eigenvalues, eigenvectors = np.linalg.eigh(hamiltonian)
    middle_pair = np.argsort(np.abs(eigenvalues))[:2]

    pair_vectors = eigenvectors[:, middle_pair]
    pair_density = np.sum(
        np.abs(pair_vectors.reshape(num_cells, 2, 2)) ** 2,
        axis=(1, 2),
    )
    pair_density /= np.sum(pair_density)

    gamma = np.diag(np.tile([1.0, -1.0], num_cells))
    chiral_error = float(
        np.linalg.norm(gamma @ hamiltonian @ gamma + hamiltonian, ord="fro")
    )
    pairing_error = float(
        np.max(np.abs(eigenvalues + eigenvalues[::-1]))
    )

    return {
        "hamiltonian": hamiltonian,
        "eigenvalues": eigenvalues,
        "eigenvectors": eigenvectors,
        "middle_pair": middle_pair,
        "pair_density": pair_density,
        "chiral_error": chiral_error,
        "pairing_error": pairing_error,
    }


def plot_results(mapping, cases):
    """Plot only the three results needed for the clean-model baseline."""
    figure, axes = plt.subplots(1, 3, figsize=(12.0, 3.8))

    # 1. TMM points and analytical SSH bands in the same coordinates.
    dispersion_axis = axes[0]
    reference = cases[0]
    q_grid = np.linspace(-np.pi, np.pi, NUM_Q)
    lower_band, upper_band = ssh_bands(
        q_grid, reference["v"], reference["w"]
    )
    dispersion_axis.plot(q_grid, upper_band, color="black", label="SSH")
    dispersion_axis.plot(q_grid, lower_band, color="black")

    allowed = mapping["allowed"]
    q_from_tmm = mapping["q"][allowed]
    lambda_from_tmm = mapping["lambda"][allowed]
    dispersion_axis.scatter(
        np.concatenate([q_from_tmm, -q_from_tmm]),
        np.tile(lambda_from_tmm, 2),
        s=5,
        color="tab:orange",
        label="TMM",
    )
    dispersion_axis.set_xlabel(r"$q=\Omega T$")
    dispersion_axis.set_ylabel(r"$\lambda$")
    dispersion_axis.set_xticks([-np.pi, 0.0, np.pi])
    dispersion_axis.set_xticklabels([r"$-\pi$", "0", r"$\pi$"])
    dispersion_axis.set_title("TMM--SSH dispersion")
    dispersion_axis.legend(fontsize=8)

    # 2. Winding paths for the two index orderings.
    winding_axis = axes[1]
    colors = ["tab:blue", "tab:green"]
    for case, color in zip(cases, colors):
        winding = int(round(case["winding"]))
        winding_axis.plot(
            case["path"].real,
            case["path"].imag,
            color=color,
            label=rf"$n_1={case['n1']:g},n_2={case['n2']:g}$: $\nu={winding}$",
        )
    winding_axis.scatter([0.0], [0.0], color="red", marker="x", s=35)
    winding_axis.axhline(0.0, color="0.75", linewidth=0.7)
    winding_axis.axvline(0.0, color="0.75", linewidth=0.7)
    winding_axis.set_aspect("equal", adjustable="box")
    winding_axis.set_xlabel(r"$\mathrm{Re}\,z$")
    winding_axis.set_ylabel(r"$\mathrm{Im}\,z$")
    winding_axis.set_title("Winding number")
    winding_axis.legend(fontsize=7)

    # 3. Open-boundary spectra for the same two cases.
    spectrum_axis = axes[2]
    for case, color in zip(cases, colors):
        eigenvalues = case["chain"]["eigenvalues"]
        spectrum_axis.scatter(
            np.arange(eigenvalues.size),
            eigenvalues,
            s=13,
            facecolors="none",
            edgecolors=color,
            label=rf"$n_1={case['n1']:g},n_2={case['n2']:g}$",
        )
    spectrum_axis.axhline(0.0, color="0.65", linewidth=0.7)
    spectrum_axis.set_xlabel("Eigenvalue index")
    spectrum_axis.set_ylabel(r"$\lambda$")
    spectrum_axis.set_title("Open-boundary spectrum")
    spectrum_axis.legend(fontsize=7)

    for axis in axes:
        axis.grid(alpha=0.2)

    figure.tight_layout()
    figure.savefig(OUTPUT_FIGURE, dpi=180)
    plt.close(figure)


def run_regression_checks(mapping, cases):
    """Fail loudly if any clean-limit identity or topology label is wrong."""
    v, w = ssh_couplings(REFERENCE_N1, REFERENCE_N2)
    q_grid = np.linspace(-np.pi, np.pi, 301)
    analytical_lower, analytical_upper = ssh_bands(q_grid, v, w)
    numerical_bands = np.array(
        [np.linalg.eigvalsh(ssh_bloch_hamiltonian(q, v, w)) for q in q_grid]
    )

    metrics = {
        "tmm_ssh_identity_error": float(
            np.max(
                np.abs(
                    mapping["trace_half"] - mapping["trace_half_from_ssh"]
                )
            )
        ),
        "ssh_bloch_band_error": float(
            np.max(
                np.abs(
                    numerical_bands
                    - np.column_stack([analytical_lower, analytical_upper])
                )
            )
        ),
    }

    expected_windings = [1.0, 0.0]
    for case, expected in zip(cases, expected_windings):
        label = "topological" if expected == 1.0 else "trivial"
        metrics[f"{label}_winding_error"] = abs(case["winding"] - expected)
        metrics[f"{label}_chiral_error"] = case["chain"]["chiral_error"]
        metrics[f"{label}_spectral_pairing_error"] = case["chain"][
            "pairing_error"
        ]
        metrics[f"{label}_minimum_abs_eigenvalue"] = float(
            np.min(np.abs(case["chain"]["eigenvalues"]))
        )

    if metrics["tmm_ssh_identity_error"] > 1e-12:
        raise AssertionError("TMM and SSH trace identities do not agree.")
    if metrics["ssh_bloch_band_error"] > 1e-12:
        raise AssertionError("Numerical and analytical SSH bands do not agree.")
    if metrics["topological_winding_error"] > 1e-12:
        raise AssertionError("n1 > n2 did not produce winding number 1.")
    if metrics["trivial_winding_error"] > 1e-12:
        raise AssertionError("n1 < n2 did not produce winding number 0.")
    if metrics["topological_chiral_error"] > 1e-12:
        raise AssertionError("The topological open chain lost chiral symmetry.")
    if metrics["trivial_chiral_error"] > 1e-12:
        raise AssertionError("The trivial open chain lost chiral symmetry.")
    if metrics["topological_minimum_abs_eigenvalue"] > 1e-8:
        raise AssertionError("The topological chain has no near-zero edge pair.")

    trivial_gap = abs(cases[1]["v"] - cases[1]["w"])
    if metrics["trivial_minimum_abs_eigenvalue"] <= trivial_gap:
        raise AssertionError("The trivial chain unexpectedly has an in-gap state.")

    return metrics


def main():
    mapping = calculate_tmm_mapping(REFERENCE_N1, REFERENCE_N2, T_BAR)

    parameter_pairs = [
        (REFERENCE_N1, REFERENCE_N2),
        (REFERENCE_N2, REFERENCE_N1),
    ]
    cases = []
    for n1, n2 in parameter_pairs:
        v, w = ssh_couplings(n1, n2)
        _, path, winding = winding_number(v, w)
        cases.append(
            {
                "n1": n1,
                "n2": n2,
                "v": v,
                "w": w,
                "path": path,
                "winding": winding,
                "chain": open_chain_data(NUM_CELLS, v, w),
            }
        )

    metrics = run_regression_checks(mapping, cases)
    plot_results(mapping, cases)

    print("Clean PTC--SSH checks passed.")
    for name, value in metrics.items():
        print(f"{name}: {value:.6e}")
    print(f"Figure: {OUTPUT_FIGURE.resolve()}")


if __name__ == "__main__":
    main()
