import apps_header

from geometry.Circle import Circle

if __name__ == "__main__":

    circle = Circle(radius=5, center=(0, 0))
    circle.set_panels(num_panels=20)
    circle.print_info()
    circle.set_output_directory("output/source_circle")
    circle.plot_circle()

    circle.set_panel_lengths()
    circle.set_panel_normals()
    circle.set_panel_tangents()

    circle.plot_normals(scale=1.0)
    circle.plot_tangents(scale=1.0)

    circle.set_control_points()
    circle.plot_control_points()
    circle.plot_panels(control_points=True)