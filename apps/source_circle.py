import apps_header

from geometry.Circle import Circle
from geometry.Shape import Shape

if __name__ == "__main__":

    circle = Circle(radius=5, center=(0, 0))
    circle.set_panels(num_panels=12)
    circle.print_info()
    circle.plot_circle()