"""Validate the exact mapping from a disordered chiral PTC to SSH.

For a short open chain, this script constructs

    K D = lambda S D

and

    H_eff psi = lambda psi,
    H_eff = S^{-1/2} K S^{-1/2},
    psi = S^{1/2} D.

It compares both spectra, checks chiral symmetry and spectral pairing, and
restores one physical electric-displacement eigenstate from ``psi``.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eigh


num_cells = 8
random_seed = 20260807
t_bar = 1e-15

# A modest disorder level is enough for this algebraic validation.  Both
# layer types are varied here because the exact mapping allows arbitrary
# positive refractive-index sequences.
n_1_0 = 2.0
n_2_0 = 2.1
disorder_fraction = 0.18
local_perturbation_fraction = 0.20

output_figure = Path(__file__).with_name(
    "disordered_ptc_ssh_mapping.png"
)
output_data = Path(__file__).with_name(
    "disordered_ptc_ssh_mapping.npz"
)


def generate_disordered_layers(
    num_cells=num_cells,
    random_seed=random_seed,
):
    """Return one reproducible positive-index chain.

    The open chain contains N internal layers ``n_1`` and N+1 external
    layers ``n_2``.  The first and last entries of ``n_2`` are the left and
    right boundary layers; ``n_2[1:-1]`` are the intercell layers.
    """
    if not isinstance(num_cells, int) or num_cells < 2:
        raise ValueError("num_cells must be an integer of at least 2.")

    rng = np.random.default_rng(random_seed)
    n_1 = n_1_0 * (
        1.0 + disorder_fraction * rng.uniform(-1.0, 1.0, num_cells)
    )
    n_2 = n_2_0 * (
        1.0
        + disorder_fraction * rng.uniform(-1.0, 1.0, num_cells + 1)
    )

    if np.any(n_1 <= 0.0) or np.any(n_2 <= 0.0):
        raise ValueError("All refractive indices must be positive.")

    # The exact chiral mapping requires tau_j / n_j = t_bar in every layer.
    tau_1 = t_bar * n_1
    tau_2 = t_bar * n_2
    return n_1, n_2, tau_1, tau_2


def build_generalized_problem(n_1, n_2):
    """Construct K, S and H_eff in the order A_1, B_1, ..., A_N, B_N."""
    n_1 = np.asarray(n_1, dtype=float)
    n_2 = np.asarray(n_2, dtype=float)

    if n_1.ndim != 1 or n_2.ndim != 1:
        raise ValueError("n_1 and n_2 must be one-dimensional arrays.")
    if n_2.size != n_1.size + 1:
        raise ValueError("N cells require N values in n_1 and N+1 in n_2.")
    if not np.all(np.isfinite(n_1)) or not np.all(np.isfinite(n_2)):
        raise ValueError("All refractive indices must be finite.")
    if np.any(n_1 <= 0.0) or np.any(n_2 <= 0.0):
        raise ValueError("All refractive indices must be positive.")

    num_cells = n_1.size
    num_sites = 2 * num_cells

    # Chronological layer order:
    # n_2,left, n_1,1, n_2,1, ..., n_1,N, n_2,right.
    n_j = np.empty(2 * num_cells + 1, dtype=float)
    n_j[0::2] = n_2
    n_j[1::2] = n_1
    r_j = 1.0 / n_j

    # Interface-site weights s_j = r_{j-1} + r_j.
    s_j = r_j[:-1] + r_j[1:]
    S = np.diag(s_j)

    # The inner layers r_j[1:-1] connect adjacent interface sites.
    K = np.zeros((num_sites, num_sites), dtype=float)
    site = np.arange(num_sites - 1)
    K[site, site + 1] = r_j[1:-1]
    K[site + 1, site] = r_j[1:-1]

    inverse_sqrt_s = 1.0 / np.sqrt(s_j)
    H_eff = (
        inverse_sqrt_s[:, np.newaxis]
        * K
        * inverse_sqrt_s[np.newaxis, :]
    )

    # Explicit unit-cell quantities from PROJECT.md.
    r_1 = 1.0 / n_1
    r_2 = 1.0 / n_2
    s_A = s_j[0::2]
    s_B = s_j[1::2]
    v = r_1 / np.sqrt(s_A * s_B)
    w = r_2[1:-1] / np.sqrt(s_B[:-1] * s_A[1:])

    return {
        "n_j": n_j,
        "r_j": r_j,
        "s_j": s_j,
        "r_1": r_1,
        "r_2": r_2,
        "s_A": s_A,
        "s_B": s_B,
        "v": v,
        "w": w,
        "K": K,
        "S": S,
        "H_eff": H_eff,
    }


def solve_and_check(mapping):
    """Solve both eigenproblems and return the required error metrics."""
    K = mapping["K"]
    S = mapping["S"]
    H_eff = mapping["H_eff"]
    s_j = mapping["s_j"]
    num_sites = H_eff.shape[0]

    # scipy.linalg.eigh handles the symmetric generalized problem directly.
    eigenvalues_generalized, D_generalized = eigh(K, S)
    eigenvalues_effective, psi = np.linalg.eigh(H_eff)

    Gamma = np.diag(np.where(np.arange(num_sites) % 2 == 0, 1.0, -1.0))

    eigenvalue_error = float(
        np.max(np.abs(eigenvalues_generalized - eigenvalues_effective))
    )
    chiral_error = float(
        np.linalg.norm(Gamma @ H_eff @ Gamma + H_eff, ord="fro")
    )
    pairing_error = float(
        np.max(np.abs(eigenvalues_effective + eigenvalues_effective[::-1]))
    )
    symmetry_error = float(
        np.linalg.norm(H_eff - H_eff.T, ord="fro")
    )
    zero_diagonal_error = float(np.max(np.abs(np.diag(H_eff))))

    row, column = np.indices(H_eff.shape)
    long_range_mask = np.abs(row - column) > 1
    long_range_error = float(np.max(np.abs(H_eff[long_range_mask])))

    # Check the explicit v_n and w_n formulas against the matrix entries.
    A_sites = np.arange(0, num_sites, 2)
    B_sites = np.arange(1, num_sites, 2)
    coupling_error = float(
        max(
            np.max(np.abs(H_eff[A_sites, B_sites] - mapping["v"])),
            np.max(
                np.abs(
                    H_eff[B_sites[:-1], A_sites[1:]] - mapping["w"]
                )
            ),
        )
    )

    # Restore the eigenstate closest to lambda = 0.
    selected_mode = int(np.argmin(np.abs(eigenvalues_effective)))
    selected_eigenvalue = float(eigenvalues_effective[selected_mode])
    psi_mode = psi[:, selected_mode]
    D_mode = psi_mode / np.sqrt(s_j)

    generalized_residual = float(
        np.linalg.norm(K @ D_mode - selected_eigenvalue * S @ D_mode)
    )

    # Generalized eigenvectors are S-normalized and have an arbitrary sign.
    D_reference = D_generalized[:, selected_mode].copy()
    if np.dot(D_reference, S @ D_mode) < 0.0:
        D_reference *= -1.0
    field_reconstruction_error = float(
        np.linalg.norm(D_reference - D_mode)
    )

    metrics = {
        "symmetry_error": symmetry_error,
        "zero_diagonal_error": zero_diagonal_error,
        "long_range_error": long_range_error,
        "coupling_error": coupling_error,
        "chiral_error": chiral_error,
        "pairing_error": pairing_error,
        "eigenvalue_error": eigenvalue_error,
        "generalized_residual": generalized_residual,
        "field_reconstruction_error": field_reconstruction_error,
    }

    return {
        "eigenvalues_generalized": eigenvalues_generalized,
        "eigenvalues_effective": eigenvalues_effective,
        "selected_mode": selected_mode,
        "selected_eigenvalue": selected_eigenvalue,
        "psi_mode": psi_mode,
        "D_mode": D_mode,
        "metrics": metrics,
    }


def calculate_local_perturbation(num_cells=num_cells):
    """Show how one external layer changes correlated effective couplings."""
    n_1_clean = np.full(num_cells, n_1_0, dtype=float)
    n_2_clean = np.full(num_cells + 1, n_2_0, dtype=float)

    # n_2[m] lies between B_m and A_{m+1}; n_2[0] is the left boundary.
    perturbed_external_layer = num_cells // 2
    n_2_perturbed = n_2_clean.copy()
    n_2_perturbed[perturbed_external_layer] *= (
        1.0 + local_perturbation_fraction
    )

    clean_mapping = build_generalized_problem(n_1_clean, n_2_clean)
    perturbed_mapping = build_generalized_problem(
        n_1_clean,
        n_2_perturbed,
    )

    # Bond order: v_1, w_1, v_2, w_2, ..., w_{N-1}, v_N.
    clean_couplings = np.empty(2 * num_cells - 1, dtype=float)
    perturbed_couplings = np.empty(2 * num_cells - 1, dtype=float)
    clean_couplings[0::2] = clean_mapping["v"]
    clean_couplings[1::2] = clean_mapping["w"]
    perturbed_couplings[0::2] = perturbed_mapping["v"]
    perturbed_couplings[1::2] = perturbed_mapping["w"]

    bond_labels = []
    for cell in range(1, num_cells + 1):
        bond_labels.append(rf"$v_{{{cell}}}$")
        if cell < num_cells:
            bond_labels.append(rf"$w_{{{cell}}}$")

    return {
        "n_1_clean": n_1_clean,
        "n_2_clean": n_2_clean,
        "n_2_perturbed": n_2_perturbed,
        "perturbed_external_layer": perturbed_external_layer,
        "bond_labels": bond_labels,
        "clean_couplings": clean_couplings,
        "perturbed_couplings": perturbed_couplings,
    }


def plot_results(mapping, results, perturbation):
    """Plot the matrix, spectra, and local correlated coupling response."""
    H_eff = mapping["H_eff"]
    num_sites = H_eff.shape[0]
    eigenvalue_index = np.arange(1, num_sites + 1)

    figure, axes = plt.subplots(1, 3, figsize=(12.0, 3.8))

    color_limit = np.max(np.abs(H_eff))
    image = axes[0].imshow(
        H_eff,
        origin="lower",
        cmap="RdBu_r",
        vmin=-color_limit,
        vmax=color_limit,
    )
    axes[0].set_xlabel("Interface site")
    axes[0].set_ylabel("Interface site")
    axes[0].set_title(r"$H_{\mathrm{eff}}$")
    figure.colorbar(image, ax=axes[0], fraction=0.046, pad=0.04)

    axes[1].scatter(
        eigenvalue_index,
        results["eigenvalues_generalized"],
        s=45,
        facecolors="none",
        edgecolors="tab:blue",
        label=r"$K D=\lambda S D$",
    )
    axes[1].scatter(
        eigenvalue_index,
        results["eigenvalues_effective"],
        s=20,
        color="tab:orange",
        marker="x",
        label=r"$H_{\mathrm{eff}}\psi=\lambda\psi$",
    )
    axes[1].axhline(0.0, color="0.6", linewidth=0.8)
    axes[1].set_xlabel("Eigenvalue index")
    axes[1].set_ylabel(r"$\lambda$")
    axes[1].set_title("Generalized--standard spectra")
    axes[1].legend(fontsize=8)

    bond_index = np.arange(1, 2 * num_cells)
    axes[2].plot(
        bond_index,
        perturbation["clean_couplings"],
        marker="o",
        linestyle="--",
        color="0.45",
        label="Clean",
    )
    axes[2].plot(
        bond_index,
        perturbation["perturbed_couplings"],
        marker="s",
        color="tab:red",
        label="One perturbed layer",
    )
    perturbed_external_layer = perturbation["perturbed_external_layer"]
    axes[2].axvline(
        2 * perturbed_external_layer,
        color="tab:blue",
        linewidth=1.0,
        linestyle=":",
        label=rf"perturbed $n_{{2,{perturbed_external_layer}}}$",
    )
    axes[2].set_xticks(bond_index)
    axes[2].set_xticklabels(
        perturbation["bond_labels"],
        rotation=60,
        fontsize=7,
    )
    axes[2].set_xlabel("Effective bond")
    axes[2].set_ylabel("Coupling")
    axes[2].set_title("Local index change: correlated bonds")
    axes[2].legend(fontsize=8)

    for axis in axes:
        axis.grid(alpha=0.2)

    figure.tight_layout()
    figure.savefig(output_figure, dpi=180)
    plt.close(figure)


def run_regression_checks(results):
    """Fail if any required identity is outside numerical precision."""
    tolerances = {
        "symmetry_error": 1e-13,
        "zero_diagonal_error": 1e-14,
        "long_range_error": 1e-14,
        "coupling_error": 1e-14,
        "chiral_error": 1e-13,
        "pairing_error": 1e-12,
        "eigenvalue_error": 1e-12,
        "generalized_residual": 1e-12,
        "field_reconstruction_error": 1e-11,
    }

    for name, value in results["metrics"].items():
        print(f"{name}: {value:.3e}")
        if value > tolerances[name]:
            raise AssertionError(
                f"{name}={value:.3e} exceeds {tolerances[name]:.1e}."
            )


def main():
    n_1, n_2, tau_1, tau_2 = generate_disordered_layers()
    mapping = build_generalized_problem(n_1, n_2)
    results = solve_and_check(mapping)
    perturbation = calculate_local_perturbation()

    run_regression_checks(results)
    plot_results(mapping, results, perturbation)

    np.savez(
        output_data,
        random_seed=random_seed,
        t_bar=t_bar,
        n_1=n_1,
        n_2=n_2,
        tau_1=tau_1,
        tau_2=tau_2,
        K=mapping["K"],
        S=mapping["S"],
        H_eff=mapping["H_eff"],
        s_j=mapping["s_j"],
        v=mapping["v"],
        w=mapping["w"],
        eigenvalues_generalized=results["eigenvalues_generalized"],
        eigenvalues_effective=results["eigenvalues_effective"],
        psi_mode=results["psi_mode"],
        D_mode=results["D_mode"],
        local_perturbation_fraction=local_perturbation_fraction,
        n_1_clean=perturbation["n_1_clean"],
        n_2_clean=perturbation["n_2_clean"],
        n_2_perturbed=perturbation["n_2_perturbed"],
        clean_couplings=perturbation["clean_couplings"],
        perturbed_couplings=perturbation["perturbed_couplings"],
    )

    print("All exact disordered-mapping checks passed.")
    print(f"Figure: {output_figure.resolve()}")
    print(f"Data: {output_data.resolve()}")




if __name__ == "__main__":
    main()
