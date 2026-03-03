"""
05b_som_modelled_response.py
Visualise the SOM response curves as implemented in the Living Dunes model.

Reads the gaussian_boost parameters directly from the ammophila_arenaria species
configuration and plots the response multiplier across the SOM gradient (0-4% SOM),
annotating with the experimental TOC range (Van Bemmelen conversion).

Output
------
    D:/projects/susana/output/som_response/som_response_curves.png
    D:/projects/susana/output/som_response/som_response_curves.pdf

Usage
-----
    python scripts/05b_som_modelled_response.py [--output-dir <path>]

This script is part of the SUSANA analysis pipeline (run_pipeline.R).
Figure copying to the manuscript directory is handled by 06_copy_figures_to_manuscript.R.
"""

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# ---------------------------------------------------------------------------
# Configuration paths
# ---------------------------------------------------------------------------
SPECIES_CONFIG_PATH = Path(
    r"C:\Users\frevdael\Documents\02_CODE\01_github\coastal_research"
    r"\living_dunes\config\default\atlantic_europe\species"
    r"\ammophila_arenaria.json"
)

DEFAULT_OUTPUT_DIR = Path(r"D:\projects\susana\output\som_response")

VAN_BEMMELEN_FACTOR: float = 1.724

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def gaussian_boost(
    x: np.ndarray,
    mean: float,
    std: float,
    amplitude: float,
    reference: float | None = None,
) -> np.ndarray:
    """Reproduce the anchored gaussian_boost function from response_functions.py.

    When reference is provided, the curve is divided by its own value at the
    reference point so that f(reference) = 1.0 exactly.
    """
    raw = 1.0 + amplitude * np.exp(-((x - mean) ** 2) / (2 * std ** 2))
    if reference is not None:
        ref_val = 1.0 + amplitude * np.exp(-((reference - mean) ** 2) / (2 * std ** 2))
        return raw / ref_val
    return raw


def load_som_factors(config_path: Path) -> dict:
    """Extract all SOM response factors from the species configuration.

    Returns a dict keyed by demographic process name, each value being
    the response_function parameters dict plus metadata.
    """
    with open(config_path, "r", encoding="utf-8") as fh:
        config = json.load(fh)

    # Navigate to the adult life phase parameters.
    # The species config wraps everything under species[0].life_phase_parameters.adult
    species_list = config.get("species", [])
    if not species_list:
        raise ValueError("No species found in config")
    adult = species_list[0]["life_phase_parameters"]["adult"]

    factors: dict = {}

    # Walk the demographic processes that may carry SOM factors
    process_paths = [
        ("height_growth",           ["environmental_factors"]),
        ("lateral_spread.kernel",   ["kernel", "environmental_factors"]),
        ("lateral_spread.target",   ["target", "environmental_factors"]),
        ("shoot_formation.kernel",  ["kernel", "environmental_factors"]),
        ("shoot_formation.target",  ["target", "environmental_factors"]),
    ]

    for label, keys in process_paths:
        process_name = label.split(".")[0]
        node = adult.get(process_name, {})
        for key in keys:
            node = node.get(key, {})
        som_cfg = node.get("soil_organic_matter", {})
        if som_cfg and "response_function" in som_cfg:
            factors[label] = {
                "enabled": som_cfg.get("enabled", False),
                "function": som_cfg["response_function"]["function"],
                "parameters": som_cfg["response_function"]["parameters"],
                "description": som_cfg.get("description", ""),
            }

    return factors


def plot_response_curves(factors: dict, output_dir: Path) -> None:
    """Generate a publication-quality panel of SOM response curves."""
    som_range = np.linspace(0.0, 4.5, 500)

    # Experimental TOC range from SUSANA (0.06 - 2.16 % TOC)
    toc_exp_min = 0.061
    toc_exp_max = 2.161
    som_exp_min = toc_exp_min * VAN_BEMMELEN_FACTOR
    som_exp_max = toc_exp_max * VAN_BEMMELEN_FACTOR

    # Reference SOM value: beach sand (R) = 0.06% TOC -> 0.1051% SOM
    som_reference = 0.061 * VAN_BEMMELEN_FACTOR

    # ------------------------------------------------------------------
    # Deduplicate: lateral_spread and shoot_formation (kernel + target)
    # all share identical parameters — only 2 distinct curves exist.
    # Show each unique parameter set once; group by colour.
    # ------------------------------------------------------------------
    unique_curves: list[dict] = []
    seen_params: list[tuple] = []

    display_groups = {
        "Height growth": ("height_growth", "#2ca02c"),
        "Lateral spread (kernel)": ("lateral_spread.kernel", "#ff7f0e"),
        "Lateral spread (target)": ("lateral_spread.target", "#ff7f0e"),
        "Shoot formation (kernel)": ("shoot_formation.kernel", "#1f77b4"),
        "Shoot formation (target)": ("shoot_formation.target", "#1f77b4"),
    }

    for label, (key, colour) in display_groups.items():
        if key not in factors:
            continue
        p = factors[key]["parameters"]
        sig = (round(p["mean"], 4), round(p["std"], 4), round(p["amplitude"], 4))
        if sig not in seen_params:
            seen_params.append(sig)
            # Build combined label for identical kernel+target pairs
            base = label.replace(" (kernel)", "").replace(" (target)", "")
            unique_curves.append({
                "label": base,
                "colour": colour,
                "params": p,
                "key": key,
            })

    # ------------------------------------------------------------------
    # Compute curves
    # ------------------------------------------------------------------
    for curve in unique_curves:
        p = curve["params"]
        cfg_ref = p.get("reference")  # from species config, e.g. 0.1051
        curve["y"] = gaussian_boost(som_range, p["mean"], p["std"], p["amplitude"], reference=cfg_ref)
        curve["f_at_ref"] = gaussian_boost(
            np.array([som_reference]), p["mean"], p["std"], p["amplitude"], reference=cfg_ref
        )[0]
        curve["f_at_zero"] = gaussian_boost(
            np.array([0.0]), p["mean"], p["std"], p["amplitude"], reference=cfg_ref
        )[0]

    # ------------------------------------------------------------------
    # Figure layout — use constrained_layout to prevent axis clashes
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(
        nrows=2, ncols=1, figsize=(9, 8),
        constrained_layout=True,
    )

    def add_secondary_toc_axis(ax: plt.Axes) -> None:
        """Add a secondary x-axis showing TOC (%) on top."""
        ax2 = ax.secondary_xaxis(
            "top",
            functions=(
                lambda s: s / VAN_BEMMELEN_FACTOR,
                lambda t: t * VAN_BEMMELEN_FACTOR,
            ),
        )
        ax2.set_xlabel("Total organic carbon (% TOC)", labelpad=6)
        ax2.tick_params(labelsize=8)

    # --- Panel (a): Individual response curves ---
    ax_a = axes[0]
    for curve in unique_curves:
        ax_a.plot(
            som_range, curve["y"],
            label=f"{curve['label']}  "
                  f"(\u03bc={curve['params']['mean']:.2f}, "
                  f"\u03c3={curve['params']['std']:.2f}, "
                  f"a={curve['params']['amplitude']:.2f})",
            color=curve["colour"],
            linewidth=2.0,
        )
        # Mark f(0) with a small dot to make the x=0 intercept explicit
        ax_a.plot(0, curve["f_at_zero"], "o", color=curve["colour"],
                  markersize=4, zorder=5)

    # Experimental range shading
    ax_a.axvspan(
        som_exp_min, som_exp_max, alpha=0.08, color="grey",
        label=f"Experimental range ({toc_exp_min:.2f}\u2013{toc_exp_max:.2f}% TOC)",
    )
    ax_a.axhline(y=1.0, color="grey", linestyle="--", linewidth=0.8, alpha=0.6,
                 label="Baseline (f = 1)")
    ax_a.axvline(x=som_reference, color="steelblue", linestyle=":", linewidth=0.9,
                 alpha=0.7, label=f"Reference sand (R, {som_reference:.2f}% SOM)")

    ax_a.set_xlabel("Soil organic matter (% SOM)", labelpad=4)
    ax_a.set_ylabel("Growth multiplier")
    ax_a.text(
        0.01, 0.97, "(a)",
        transform=ax_a.transAxes, fontsize=10, fontweight="bold", va="top"
    )
    ax_a.legend(fontsize=7.5, loc="upper right", framealpha=0.9,
                handlelength=1.6, borderpad=0.6)
    ax_a.set_xlim(0, 4.5)
    ax_a.set_ylim(0.88, 2.2)
    add_secondary_toc_axis(ax_a)

    # --- Panel (b): Combined geometric mean ---
    ax_b = axes[1]
    all_y = np.stack([c["y"] for c in unique_curves], axis=0)
    geo_mean = np.prod(all_y, axis=0) ** (1.0 / len(unique_curves))

    ax_b.plot(
        som_range, geo_mean, color="#333333", linewidth=2.5,
        label=f"Geometric mean of {len(unique_curves)} SOM factors",
    )
    ax_b.fill_between(
        som_range, 1.0, geo_mean,
        where=geo_mean > 1.0, alpha=0.15, color="#2ca02c",
        label="Boost region",
    )
    ax_b.fill_between(
        som_range, geo_mean, 1.0,
        where=geo_mean < 1.0, alpha=0.15, color="#d62728",
        label="Penalty region (below baseline)",
    )
    ax_b.axvspan(som_exp_min, som_exp_max, alpha=0.08, color="grey",
                 label="Experimental range")
    ax_b.axhline(y=1.0, color="grey", linestyle="--", linewidth=0.8, alpha=0.6)
    ax_b.axvline(x=som_reference, color="steelblue", linestyle=":", linewidth=0.9, alpha=0.7)

    ax_b.set_xlabel("Soil organic matter (% SOM)", labelpad=4)
    ax_b.set_ylabel("Combined growth multiplier")
    ax_b.text(
        0.01, 0.97, "(b)",
        transform=ax_b.transAxes, fontsize=10, fontweight="bold", va="top"
    )
    ax_b.legend(fontsize=7.5, loc="upper right", framealpha=0.9)
    ax_b.set_xlim(0, 4.5)
    y_max = max(1.85, float(np.max(geo_mean)) * 1.06)
    ax_b.set_ylim(0.88, y_max)
    add_secondary_toc_axis(ax_b)

    # Notes below both panels
    note_text = (
        "Density growth SOM factor excluded (outlier-driven, p\u00a0=\u00a00.09 after removal). "
        "Curves use anchored boost: f(SOM)/f(reference) so that f(reference sand) = 1.0 exactly. "
        "Shaded band = SUSANA experimental range (0.06\u20132.16\u00a0% TOC). "
        "Blue dotted line = reference beach sand (R, 0.10\u00a0% SOM)."
    )
    fig.text(
        0.08, -0.01, note_text,
        fontsize=7, fontstyle="italic", color="grey", va="top", wrap=True,
    )

    # ------------------------------------------------------------------
    # Save — primary location: D:\projects\susana\output\som_response
    # Manuscript figure copying is handled by 06_copy_figures_to_manuscript.R
    # ------------------------------------------------------------------
    output_dir.mkdir(parents=True, exist_ok=True)

    stem = "som_response_curves"
    fig.savefig(output_dir / f"{stem}.png", dpi=300, bbox_inches="tight",
                facecolor="white")
    fig.savefig(output_dir / f"{stem}.pdf", bbox_inches="tight",
                facecolor="white")
    print(f"Saved: {output_dir / stem}.png")
    print(f"Saved: {output_dir / stem}.pdf")
    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    """Entry point."""
    parser = argparse.ArgumentParser(
        description="Visualise Living Dunes SOM response curves."
    )
    parser.add_argument(
        "--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR,
        help="Directory for output figures (default: output/som_response)."
    )
    parser.add_argument(
        "--config", type=Path, default=SPECIES_CONFIG_PATH,
        help="Path to ammophila_arenaria.json species config."
    )
    args = parser.parse_args()

    print(f"Loading species config: {args.config}")
    factors = load_som_factors(args.config)

    if not factors:
        print("No SOM factors found in species config.", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(factors)} SOM response factors:")
    for key, info in factors.items():
        status = "ENABLED" if info["enabled"] else "disabled"
        params = info["parameters"]
        print(f"  {key}: {info['function']}  mu={params['mean']:.3f}  "
              f"sigma={params['std']:.3f}  a={params['amplitude']:.3f}  [{status}]")

    plot_response_curves(factors, args.output_dir)


if __name__ == "__main__":
    main()
