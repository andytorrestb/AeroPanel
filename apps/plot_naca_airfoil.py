import apps_header

from geometry.NACA import NACA

if __name__ == "__main__":
    # Example: NACA 0012, 80 panels with cosine spacing
    airfoil = NACA("0012", chord_length=1.0)
    airfoil.set_panels(num_panels=80, spacing="cosine", angle_of_attack_deg=0.0)
    airfoil.set_output_directory("output/source_naca")

    airfoil.print_info()

    # Plot outline and diagnostics similar to Circle plots
    airfoil.plot_airfoil()
    airfoil.plot_normals(scale=0.2)
    airfoil.plot_tangents(scale=0.2)
    airfoil.plot_control_points()
