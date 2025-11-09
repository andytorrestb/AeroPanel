import os
import numpy as np
import math

import matplotlib.pyplot as plt

class Shape:
    def print_info(self):
        raise NotImplementedError("Subclasses must implement this method")
    
    def set_output_directory(self, output_dir: str):
        self.output_dir = output_dir

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    
    def set_panel_lengths(self):
        # Chord length between consecutive panel boundary points (wrap-around)
        S = np.zeros(self.num_panels)
        for i in range(self.num_panels):
            j = (i + 1) % self.num_panels
            x_curr = self.panels[i]["x_i"]
            y_curr = self.panels[i]["y_i"]
            x_next = self.panels[j]["x_i"]
            y_next = self.panels[j]["y_i"]
            S[i] = math.hypot(x_next - x_curr, y_next - y_curr)
        self.S = S

    def set_panel_normals(self):
        # For the circle implementation we adopt CW param: y = cy - r*sin(theta)
        # Outward unit normal is radial: [cos(theta), -sin(theta)]
        N = np.zeros((self.num_panels, 2))
        for i in range(self.num_panels):
            theta = math.radians(self.panels[i]["mid_angle"])
            N[i] = [math.cos(theta), -math.sin(theta)]
        self.N = N

    def set_panel_tangents(self):
        # With CW paramization, increasing angle moves clockwise.
        # Unit tangent along panel direction is d/dtheta [cos, -sin] = [-sin, -cos]
        T = np.zeros((self.num_panels, 2))
        for i in range(self.num_panels):
            theta = math.radians(self.panels[i]["mid_angle"])
            T[i] = [-math.sin(theta), -math.cos(theta)]
        self.T = T

    def plot_normals(self, scale=1.0):
        if not hasattr(self, 'N'):
            raise ValueError("Panel normals not set. Call set_panel_normals() first.")

        fig, ax = plt.subplots()
        ax.set_title("Panel Normals")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")

        for i in range(self.num_panels):
            # Prefer the panel midpoint if available
            x_i = self.panels[i].get("xm_i", self.panels[i]["x_i"]) 
            y_i = self.panels[i].get("ym_i", self.panels[i]["y_i"]) 
            nx, ny = self.N[i]
            ax.arrow(x_i, y_i, scale * nx, scale * ny,
                     head_width=0.1, head_length=0.2, fc='b', ec='b')

        ax.set_aspect('equal', 'box')
        ax.set_xlim(self.center[0] - self.radius - 1, self.center[0] + self.radius + 1)
        ax.set_ylim(self.center[1] - self.radius - 1, self.center[1] + self.radius + 1)
        ax.grid(True)
        plt.savefig(self.output_dir + "/panel_normals.png")

    def plot_tangents(self, scale=1.0):
        if not hasattr(self, 'T'):
            raise ValueError("Panel tangents not set. Call set_panel_tangents() first.")

        fig, ax = plt.subplots()
        ax.set_title("Panel Tangents")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")

        for i in range(self.num_panels):
            # Prefer the panel midpoint if available
            x_i = self.panels[i].get("xm_i", self.panels[i]["x_i"]) 
            y_i = self.panels[i].get("ym_i", self.panels[i]["y_i"]) 
            tx, ty = self.T[i]
            ax.arrow(x_i, y_i, scale * tx, scale * ty,
                     head_width=0.1, head_length=0.2, fc='g', ec='g')

        ax.set_aspect('equal', 'box')
        ax.set_xlim(self.center[0] - self.radius - 1, self.center[0] + self.radius + 1)
        ax.set_ylim(self.center[1] - self.radius - 1, self.center[1] + self.radius + 1)
        ax.grid(True)
        plt.savefig(self.output_dir + "/panel_tangents.png")


    def set_control_points(self):
        # Set control points at 3/4 of each panel length from the start point
        self.control_points = []
        for i in range(self.num_panels - 1):
            start_x = self.panels[i]["x_i"]
            end_x = self.panels[i+1]["x_i"]
            control_x = start_x + 0.75 * (end_x - start_x)

            start_y = self.panels[i]["y_i"]
            end_y = self.panels[i+1]["y_i"]
            control_y = start_y + 0.75 * (end_y - start_y)

            self.control_points.append((control_x, control_y))
        self.control_points = np.array(self.control_points)
    
    def plot_control_points(self):
        if not hasattr(self, 'control_points'):
            raise ValueError("Control points not set. Call set_control_points() first.")

        fig, ax = plt.subplots()
        ax.set_title("Panel Control Points")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")

        ax.plot(self.control_points[:, 0], self.control_points[:, 1], 'ms', label='Control Points')

        ax.set_aspect('equal', 'box')
        ax.set_xlim(self.center[0] - self.radius - 1, self.center[0] + self.radius + 1)
        ax.set_ylim(self.center[1] - self.radius - 1, self.center[1] + self.radius + 1)
        ax.grid(True)
        ax.legend()
        plt.savefig(self.output_dir + "/panel_control_points.png")
    
    def plot_panels(self, control_points=False):
        if not hasattr(self, 'panels'):
            raise ValueError("Panels not set. Call set_panels() first.")

        fig, ax = plt.subplots()
        ax.set_title("Panels")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")

        x_i = []
        y_i = []
        for panel in self.panels:
            x_start = panel["x_i"]
            y_start = panel["y_i"]
            x_i.append(x_start)
            y_i.append(y_start)
        x_i = np.array(x_i)
        y_i = np.array(y_i)
        ax.plot(x_i, y_i, 'r-')
        ax.plot(x_i, y_i, 'ro')

        if control_points:
            if not hasattr(self, 'control_points'):
                raise ValueError("Control points not set. Call set_control_points() first.")
            ax.plot(self.control_points[:, 0], self.control_points[:, 1], 'ms', label='Control Points')

        ax.set_aspect('equal', 'box')
        ax.set_xlim(self.center[0] - self.radius - 1, self.center[0] + self.radius + 1)
        ax.set_ylim(self.center[1] - self.radius - 1, self.center[1] + self.radius + 1)
        ax.grid(True)
        plt.savefig(self.output_dir + "/panels.png")