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