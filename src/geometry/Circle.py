import math

import numpy as np
import matplotlib.pyplot as plt

from geometry.Shape import Shape

class Circle(Shape):
    def __init__(self, radius: float, center: tuple = (0, 0)):
        self.radius = radius
        self.center = center
        self.num_panels = 0
        self.panel_angle = 0.0
        # Per-panel mid-point coordinates (x_i, y_i)
        self.x_i = []
        self.y_i = []
        # Detailed per-panel info (angles, mid-point)
        self.panels = []

    def set_panels(self, num_panels: int):
        """Discretize the circle into panels and compute per-panel mid-point (x_i, y_i).

        Angle convention: 0 deg is along +x axis; positive rotation is CLOCKWISE.
        """
        if num_panels <= 0:
            raise ValueError("num_panels must be a positive integer")

        self.num_panels = num_panels
        self.panel_angle = 360.0 / num_panels

        cx, cy = self.center

        # Reset containers
        self.x_i = []
        self.y_i = []
        self.panels = []

        for i in range(num_panels):
            start_deg = i * self.panel_angle
            end_deg = (i + 1) * self.panel_angle
            mid_deg = 0.5 * (start_deg + end_deg)

            # Convert to radians. To make positive rotation clockwise,
            # use y = cy - r*sin(theta) (i.e., flip the sign on sin).
            theta = math.radians(start_deg)
            xi = cx + self.radius * math.cos(theta)
            yi = cy - self.radius * math.sin(theta)

            self.x_i.append(xi)
            self.y_i.append(yi)

            self.panels.append(
                {
                    "start_angle": start_deg,
                    "end_angle": end_deg,
                    "mid_angle": mid_deg,
                    "x_i": xi,
                    "y_i": yi,
                }
            )

        self.x_i = np.array(self.x_i)
        self.y_i = np.array(self.y_i)

        # Re-append first point to close the circle when plotting
        self.x_i = np.r_[self.x_i, self.x_i[0]]
        self.y_i = np.r_[self.y_i, self.y_i[0]]
    
    def plot_circle(self):

        fig, ax = plt.subplots()

        # Plot panel start points and connect them with a continuous line
        if self.num_panels > 0:
            # Connect the points in order and close the loop to form the circle outline
            ax.plot(self.x_i, self.y_i, 'ro', label='Panel Start Points')
            ax.plot(self.x_i, self.y_i, 'r-', label='Panel Start Points')

        ax.set_aspect('equal', 'box')
        ax.set_xlim(self.center[0] - self.radius - 1, self.center[0] + self.radius + 1)
        ax.set_ylim(self.center[1] - self.radius - 1, self.center[1] + self.radius + 1)
        ax.set_title('Circle with Panels')
        ax.set_xlabel('X-axis')
        ax.set_ylabel('Y-axis')
        ax.grid(True)
        if self.num_panels > 0:
            ax.legend()
        plt.savefig("circle_with_panels.png")

    def print_info(self):
        print(f"Circle with radius: {self.radius} and center: {self.center}")
        print(f"Number of panels: {self.num_panels}")
        if not self.panels:
            print("Panels not initialized. Call set_panels(num_panels) first.")
            return

        for i, panel in enumerate(self.panels):
            print(
                f"Panel {i}: Start={panel['start_angle']:.2f}°, End={panel['end_angle']:.2f}°, "
                f"Mid={panel['mid_angle']:.2f}°, x_i={panel['x_i']:.4f}, y_i={panel['y_i']:.4f}"
            )

        
