import apps_header

import numpy as np

from geometry.Circle import Circle
from solver.PanelSourceSolver import PanelSourceSolver as pss

if __name__ == "__main__":

    N = [20, 40, 80, 160]  # Number of panels to test

    results = {}

    for n_panels in N:
        # Set up the circle geometry (+ pre-processing)
        circle = Circle(radius=5, center=(0, 0))
        circle.set_panels(num_panels=n_panels)
        circle.set_panels_cosine(num_panels=n_panels)
        circle.set_output_directory("output/source_cylinder")
        circle.set_panel_lengths()
        circle.set_panel_normals()
        circle.set_panel_tangents()
        circle.set_control_points()


        # Set up and solve the panel source problem
        solver = pss(shape=circle, U_inf=1.0)
        solver.assemble_influence_matrix()
        solver.assemble_rhs()
        solver.solve_source_strengths()
        solver.compute_tangential_velocities()
        solver.compute_pressure_coefficients()

        # Output results
        solver.print_panel_source_strengths()
        solver.print_panel_tangential_velocities()
        solver.print_panel_pressure_coefficients()

        # Plot Cp distribution
        # solver.plot_pressure_coefficients()
        theta, Cp = solver.get_Cp_theta()
        results[n_panels] = (theta, Cp)

    # Comparative Cp plot for different panel counts vs analytcal solution

    theta_analytical = np.linspace(0, 2 * np.pi, 200)
    Cp_analytical = 1 - 4 * (np.sin(theta_analytical))**2
    theta_analytical = theta_analytical * 180.0 / np.pi  # Convert to degrees for consistency


    import matplotlib.pyplot as plt
    # Cartesian comparison
    fig, ax = plt.subplots()
    ax.plot(theta_analytical, Cp_analytical, 'k--', label="Analytical Solution")
    for n_panels, (theta, Cp) in results.items():
        ax.plot(theta, Cp, label=f"{n_panels} panels")
    ax.set_xlabel("Theta (degrees)")
    ax.set_ylabel("Pressure Coefficient (Cp)")
    ax.set_title("Pressure Coefficient Distribution for Different Panel Counts")
    ax.grid(True)
    ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), borderaxespad=0.0, frameon=True)
    fig.savefig("output/source_cylinder/comparative_Cp.png", bbox_inches='tight')

    # Polar comparison (theta must be in radians)
    # Negative Cp values are projected opposite by adding π to angle and plotting |Cp| as radius.
    fig_polar, ax_polar = plt.subplots(subplot_kw={'projection': 'polar'})
    theta_analytical_rad = theta_analytical * np.pi / 180.0
    Cp_a = Cp_analytical
    theta_a_proj = np.where(Cp_a >= 0.0, theta_analytical_rad, theta_analytical_rad + np.pi)
    r_a = np.abs(Cp_a)
    ax_polar.plot(theta_a_proj, r_a, 'k--', label="Analytical Solution")
    for n_panels, (theta, Cp) in results.items():
        theta_rad = np.deg2rad(theta)
        theta_proj = np.where(Cp >= 0.0, theta_rad, theta_rad + np.pi)
        r = np.abs(Cp)
        ax_polar.plot(theta_proj, r, label=f"{n_panels} panels")
    ax_polar.set_theta_zero_location('E')  # 0 degrees at +x
    ax_polar.set_theta_direction(-1)       # CW increasing, to match geometry convention
    ax_polar.set_title("|Cp| with Sign-Projected Direction (Polar)")
    # Place legend outside to the right of the polar axes
    ax_polar.legend(loc='center left', bbox_to_anchor=(1.15, 0.5), frameon=True)
    fig_polar.savefig("output/source_cylinder/comparative_Cp_polar.png", bbox_inches='tight')