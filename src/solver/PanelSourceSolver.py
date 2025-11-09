import numpy as np

class PanelSourceSolver:
    def __init__(self, shape, U_inf):
        print
        self.shape = shape
        self.U_inf = U_inf
        return
    
    def assemble_influence_matrix(self):
        print("Assembling influence matrix...")
        self.M = np.zeros((self.shape.num_panels, self.shape.num_panels))

        for i in range(self.shape.num_panels):
            for j in range(self.shape.num_panels):
                if i == j:
                    self.M[i, j] = 0.5  # Self-influence
                else:
                    # Simplified influence calculation for demonstration
                    self.M[i, j] = 1.0 / (abs(i - j) + 1)

        print(f"Influence matrix assembled with shape {self.M.shape}.")
        return
    
    def assemble_rhs(self):
        print("Assembling RHS vector...")
        self.RHS = np.zeros(self.shape.num_panels)

        for i in range(self.shape.num_panels):
            # Simplified RHS calculation for demonstration
            self.RHS[i] = -self.U_inf * np.cos(np.radians(self.shape.panels[i]["mid_angle"]))

        print(f"RHS vector assembled with shape {self.RHS.shape}.")
        return
    
    def solve_source_strengths(self):
        print("Solving for source strengths...")
        self.sigma = np.linalg.solve(self.M, self.RHS)
        print("Source strengths computed.")
        return
    
    def compute_tangential_velocities(self):
        print("Computing tangential velocities...")
        self.V_t = np.zeros(self.shape.num_panels)

        for i in range(self.shape.num_panels):
            # Simplified tangential velocity calculation for demonstration
            self.V_t[i] = self.U_inf + np.sum(self.sigma) / (2 * np.pi * self.shape.radius)

        print("Tangential velocities computed.")
        return
    
    def compute_pressure_coefficients(self):
        print("Computing pressure coefficients...")
        self.Cp = np.zeros(self.shape.num_panels)

        for i in range(self.shape.num_panels):
            V = self.V_t[i]
            self.Cp[i] = 1 - (V / self.U_inf) ** 2

        print("Pressure coefficients computed.")
        return
    
    def print_panel_source_strengths(self):
        print("Panel Source Strengths (sigma):")
        for i in range(self.shape.num_panels):
            print(f"Panel {i + 1}: sigma = {self.sigma[i]:.4f}")
        return
    
    def print_panel_tangential_velocities(self):
        print("Panel Tangential Velocities (V_t):")
        for i in range(self.shape.num_panels):
            print(f"Panel {i + 1}: V_t = {self.V_t[i]:.4f}")
        return
    
    def print_panel_pressure_coefficients(self):
        print("Panel Pressure Coefficients (Cp):")
        for i in range(self.shape.num_panels):
            print(f"Panel {i + 1}: Cp = {self.Cp[i]:.4f}")
        return
    
    def plot_pressure_coefficients(self):
        import matplotlib.pyplot as plt

        print("Plotting pressure coefficients...")
        theta = [panel["mid_angle"] for panel in self.shape.panels]

        plt.figure()
        plt.plot(theta, self.Cp, marker='o')
        plt.xlabel("Panel Mid-Angle (degrees)")
        plt.ylabel("Pressure Coefficient (Cp)")
        plt.title("Pressure Coefficient Distribution")
        plt.gca().invert_yaxis()  # Cp is typically plotted inverted
        plt.grid(True)
        plt.savefig(self.shape.output_dir + "/pressure_coefficients.png")
        plt.close()
        print("Pressure coefficient plot saved.")
        return