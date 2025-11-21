import apps_header

from geometry.NACA import NACA

APPLY_KUTTA_CONDITION = True

if __name__ == "__main__":
    # Example: NACA 0012, 80 panels with cosine spacing
    airfoil = NACA("0012", chord_length=1.0)
    airfoil.set_output_directory("output/naca_0012")
    airfoil.set_panels(num_panels=80, spacing="cosine", angle_of_attack_deg=5.0)
    airfoil.set_panel_normals()
    airfoil.set_panel_tangents()
    airfoil.set_panel_lengths()
    airfoil.set_control_points()
    # airfoil.print_info()

    # Solver
    from solver.PanelSourceSolver import PanelSourceSolver as pss
    solver = pss(shape=airfoil, U_inf=1.0, enforce_kutta=APPLY_KUTTA_CONDITION)
    solver.assemble_influence_matrix()
    solver.assemble_rhs()
    solver.solve_source_strengths()
    solver.compute_tangential_velocities()
    solver.compute_pressure_coefficients()
    solver.compute_aero_coefficients()
    # solver.print_panel_source_strengths()
    # solver.print_panel_tangential_velocities()
    # solver.print_panel_pressure_coefficients()
    solver.plot_pressure_coefficients()
    # If you want x/c vs Cp data for further plotting:
    # x_over_c, Cp = solver.get_Cp_x()
