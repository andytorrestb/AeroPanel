import apps_header

from geometry.Circle import Circle
from solver.PanelSourceSolver import PanelSourceSolver as pss

if __name__ == "__main__":

    # Set up the circle geometry (+ pre-processing)
    circle = Circle(radius=5, center=(0, 0))
    circle.set_panels(num_panels=20)
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
    solver.plot_pressure_coefficients()