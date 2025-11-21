import csv
from pathlib import Path

import apps_header
import aerosandbox as asb
import matplotlib.pyplot as plt
import numpy as np

from geometry.NACA import NACA


APPLY_KUTTA_CONDITION = True
ANGLE_OF_ATTACK_DEG = 5.0
XFOIL_VISCOSITY = True
XFOIL_REYNOLDS = 1_000_000.0


def extract_xfoil_cp_surfaces(boundary_layer_df):
    """Split XFOIL boundary-layer dump into upper/lower Cp arrays."""
    x = boundary_layer_df["x"].to_numpy(dtype=float)
    y = boundary_layer_df["y"].to_numpy(dtype=float)
    cp = boundary_layer_df["Cp"].to_numpy(dtype=float)

    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(cp)
    x = x[mask]
    y = y[mask]
    cp = cp[mask]

    # Remove wake samples that extend beyond the trailing edge (XFoil reports x>1 in the wake)
    wake_mask = (np.abs(y) <= 1e-5) & (x > 1.01)
    if np.any(wake_mask):
        x = x[~wake_mask]
        y = y[~wake_mask]
        cp = cp[~wake_mask]

    # Normalize x so the chord matches the unit-chord panel method solution
    x_min = np.min(x)
    x_max = np.max(x)
    chord = x_max - x_min
    if chord <= 0.0:
        chord = 1.0
    x = (x - x_min) / chord
    x = np.clip(x, 0.0, 1.0)

    upper_mask = y >= 0.0
    lower_mask = ~upper_mask

    upper_sorted = np.argsort(x[upper_mask])
    lower_sorted = np.argsort(x[lower_mask])

    upper_data = {
        "x": x[upper_mask][upper_sorted],
        "Cp": cp[upper_mask][upper_sorted],
    }
    lower_data = {
        "x": x[lower_mask][lower_sorted],
        "Cp": cp[lower_mask][lower_sorted],
    }

    return {"upper": upper_data, "lower": lower_data}


def write_cp_csv(filepath, cp_data):
    """Persist Cp surface data to CSV for external plotting."""
    with open(filepath, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["surface", "x_over_c", "Cp"])
        for surface_name in ("upper", "lower"):
            surface = cp_data.get(surface_name, {})
            x_values = surface.get("x", [])
            cp_values = surface.get("Cp", [])
            for x_val, cp_val in zip(x_values, cp_values):
                writer.writerow([surface_name, f"{x_val:.6f}", f"{cp_val:.6f}"])


def plot_cp_comparison(panel_data, xfoil_data, output_dir):
    """Generate Cp distribution comparison plot between solvers."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(10, 6))
    plt.plot(
        panel_data["upper"]["x"],
        panel_data["upper"]["Cp"],
        "b-o",
        label="Panel Upper",
        linewidth=2,
        markersize=4,
    )
    plt.plot(
        panel_data["lower"]["x"],
        panel_data["lower"]["Cp"],
        "r-o",
        label="Panel Lower",
        linewidth=2,
        markersize=4,
    )
    plt.plot(
        xfoil_data["upper"]["x"],
        xfoil_data["upper"]["Cp"],
        "b--",
        label="XFOIL Upper",
        linewidth=2,
    )
    plt.plot(
        xfoil_data["lower"]["x"],
        xfoil_data["lower"]["Cp"],
        "r--",
        label="XFOIL Lower",
        linewidth=2,
    )
    plt.xlabel("x / c")
    plt.ylabel("Pressure Coefficient (Cp)")
    plt.title("Panel Method vs. XFOIL Cp Comparison")
    plt.gca().invert_yaxis()
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    output_path = output_dir / "pressure_coefficients_comparison.png"
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Comparison plot saved to {output_path}")


if __name__ == "__main__":
    # Example: NACA 0012, 80 panels with cosine spacing
    airfoil = NACA("0012", chord_length=1.0)
    airfoil.set_output_directory("output/xfoil_source_comp")
    airfoil.set_panels(
        num_panels=80,
        spacing="cosine",
        angle_of_attack_deg=-1 * ANGLE_OF_ATTACK_DEG,
    )
    airfoil.set_panel_normals()
    airfoil.set_panel_tangents()
    airfoil.set_panel_lengths()
    airfoil.set_control_points()

    # Source panel Solver
    from solver.PanelSourceSolver import PanelSourceSolver as pss

    solver = pss(shape=airfoil, U_inf=1.0, enforce_kutta=APPLY_KUTTA_CONDITION)
    solver.assemble_influence_matrix()
    solver.assemble_rhs()
    solver.solve_source_strengths()
    solver.compute_tangential_velocities()
    solver.compute_pressure_coefficients()
    solver.compute_aero_coefficients()
    solver.plot_pressure_coefficients()
    panel_cp_data = solver.get_Cp_surfaces()

    # XFOIL source panel solver for comparison
    if not XFOIL_VISCOSITY:
        raise ValueError(
            "Cp extraction requires viscous XFOIL mode; set XFOIL_VISCOSITY to True."
        )

    xfoil_airfoil = asb.Airfoil("NACA0012").repanel(n_points_per_side=200)
    xfoil_solver = asb.XFoil(
        airfoil=xfoil_airfoil,
        Re=XFOIL_REYNOLDS,
        mach=0.0,
        max_iter=200,
        include_bl_data=True,
        xfoil_repanel=True,
        xfoil_command="/mnt/c/Users/andy9/Desktop/xfoil/xfoil.exe",  # or the full Windows path if running from Windows Python
    )

    xfoil_results = xfoil_solver.alpha(ANGLE_OF_ATTACK_DEG)
    if len(xfoil_results.get("alpha", [])) == 0:
        raise RuntimeError("XFOIL did not return any converged solutions.")

    target_idx = int(np.argmin(np.abs(xfoil_results["alpha"] - ANGLE_OF_ATTACK_DEG)))
    cl_xfoil = float(xfoil_results["CL"][target_idx])
    cd_xfoil = float(xfoil_results["CD"][target_idx])
    cm_xfoil = float(xfoil_results["CM"][target_idx])

    boundary_layers = xfoil_results.get("bl_data")
    if boundary_layers is None:
        raise RuntimeError(
            "XFOIL run did not return boundary-layer data; cannot extract Cp distribution."
        )

    xfoil_cp_data = extract_xfoil_cp_surfaces(boundary_layers[target_idx])

    output_directory = Path(airfoil.output_dir)
    output_directory.mkdir(parents=True, exist_ok=True)
    panel_csv_path = output_directory / "panel_method_cp.csv"
    xfoil_csv_path = output_directory / "xfoil_cp.csv"
    write_cp_csv(panel_csv_path, panel_cp_data)
    write_cp_csv(xfoil_csv_path, xfoil_cp_data)
    print(f"Panel-method Cp saved to {panel_csv_path}")
    print(f"XFOIL Cp saved to {xfoil_csv_path}")
    plot_cp_comparison(panel_cp_data, xfoil_cp_data, output_directory)

    print("\nPanel Method Aerodynamics:")
    print(f"  Cl = {solver.Cl:.4f}")
    print(f"  Cm_c/4 = {solver.Cm_c4:.4f}")

    print("\nXFOIL Aerodynamics:")
    print(f"  alpha = {xfoil_results['alpha'][target_idx]:.3f} deg")
    print(f"  Cl = {cl_xfoil:.4f}")
    print(f"  Cd = {cd_xfoil:.5f}")
    print(f"  Cm = {cm_xfoil:.4f}")

    print("\nComparative Metrics:")
    print(f"  Cl difference (panel - XFOIL) = {solver.Cl - cl_xfoil:+.4f}")