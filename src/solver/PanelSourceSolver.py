import numpy as np

class PanelSourceSolver:
    def __init__(self, shape, U_inf):
        print
        self.shape = shape
        self.U_inf = U_inf
        return
    
    def assemble_influence_matrix(self):
        print("Assembling influence matrix...")
        Np = self.shape.num_panels
        self.M = np.zeros((Np, Np))

        # Precompute easy accessors
        normals = getattr(self.shape, 'N', None)
        if normals is None:
            raise ValueError("Panel normals not set on shape. Call shape.set_panel_normals() before assembling.")

        # For endpoints of panel i: start is panels[i], end is panels[(i+1) % Np]
        def panel_endpoints(i):
            i_next = (i + 1) % Np
            x1 = self.shape.panels[i]["x_i"]
            y1 = self.shape.panels[i]["y_i"]
            x2 = self.shape.panels[i_next]["x_i"]
            y2 = self.shape.panels[i_next]["y_i"]
            return x1, y1, x2, y2

        # Control point for row j: use panel midpoint (xm, ym)
        def control_point(j):
            return self.shape.panels[j].get("xm_i", self.shape.panels[j]["x_i"]), \
                   self.shape.panels[j].get("ym_i", self.shape.panels[j]["y_i"])

        for j in range(Np):
            xP, yP = control_point(j)
            nj = normals[j]  # outward unit normal at control point j

            for i in range(Np):
                if i == j:
                    # Analytical self-influence of a source panel is 1/2 on the normal component
                    self.M[j, i] = 0.5
                    continue

                x1, y1, x2, y2 = panel_endpoints(i)
                u_vec = self._source_panel_velocity_at_point(x1, y1, x2, y2, xP, yP)
                # Project the induced velocity onto control point normal: A_ji = n_j · u^(S) for σ_i = 1
                self.M[j, i] = np.dot(nj, u_vec)

        print(f"Influence matrix assembled with shape {self.M.shape}.")
        return

    @staticmethod
    def _source_panel_velocity_at_point(x1, y1, x2, y2, xP, yP):
        """Velocity induced at point P by a straight source panel of unit strength (σ=1).

        Uses standard analytic integrals in the local (panel-aligned) frame and rotates back
        to global coordinates. Returns [u, v].
        """
        dx = x2 - x1
        dy = y2 - y1
        S = np.hypot(dx, dy)
        if S == 0.0:
            return np.array([0.0, 0.0])

        beta = np.arctan2(dy, dx)  # panel angle relative to global x-axis
        cosb = np.cos(beta)
        sinb = np.sin(beta)

        # Transform field point into panel coordinate system (x', y')
        x_rel = xP - x1
        y_rel = yP - y1
        x_p = x_rel * cosb + y_rel * sinb
        y_p = -x_rel * sinb + y_rel * cosb

        # Distances to panel ends in local coordinates
        r1_sq = x_p**2 + y_p**2
        r2_sq = (x_p - S)**2 + y_p**2
        r1 = np.sqrt(max(r1_sq, 1e-30))
        r2 = np.sqrt(max(r2_sq, 1e-30))

        # Angles from the local x'-axis to vectors to the endpoints
        # Add small epsilon to avoid undefined atan2 when y_p=0 and x'=0
        eps = 1e-15
        theta1 = np.arctan2(y_p, x_p + eps)
        theta2 = np.arctan2(y_p, (x_p - S) + eps)

        # Induced velocity components in local frame (unit-strength source)
        u_xp = (1.0 / (2.0 * np.pi)) * np.log(r2 / r1)
        u_yp = (1.0 / (2.0 * np.pi)) * (theta2 - theta1)

        # Rotate back to global coordinates
        u = u_xp * cosb - u_yp * sinb
        v = u_xp * sinb + u_yp * cosb
        return np.array([u, v])
    
    def assemble_rhs(self):
        print("Assembling RHS vector...")
        Np = self.shape.num_panels
        self.RHS = np.zeros(Np)

        # Freestream assumed along +x direction with magnitude U_inf
        U_vec = np.array([self.U_inf, 0.0])
        normals = getattr(self.shape, 'N', None)
        if normals is None:
            raise ValueError("Panel normals not set on shape. Call shape.set_panel_normals().")

        for j in range(Np):
            n_j = normals[j]
            # RHS_j = - n_j · U_inf
            self.RHS[j] = -np.dot(n_j, U_vec)

        print(f"RHS vector assembled with shape {self.RHS.shape}.")
        return
    
    def solve_source_strengths(self):
        print("Solving for source strengths...")
        self.sigma = np.linalg.solve(self.M, self.RHS)
        print("Source strengths computed.")
        return
    
    def compute_tangential_velocities(self):
        print("Computing tangential velocities...")
        Np = self.shape.num_panels
        tangents = getattr(self.shape, 'T', None)
        normals = getattr(self.shape, 'N', None)
        if tangents is None or normals is None:
            raise ValueError("Panel tangents/normals not set. Call shape.set_panel_tangents() and set_panel_normals().")

        self.V_t = np.zeros(Np)
        U_vec = np.array([self.U_inf, 0.0])

        # Precompute panel endpoint data for reuse
        def panel_endpoints(i):
            i_next = (i + 1) % Np
            x1 = self.shape.panels[i]["x_i"]
            y1 = self.shape.panels[i]["y_i"]
            x2 = self.shape.panels[i_next]["x_i"]
            y2 = self.shape.panels[i_next]["y_i"]
            return x1, y1, x2, y2

        def control_point(j):
            return self.shape.panels[j].get("xm_i", self.shape.panels[j]["x_i"]), \
                   self.shape.panels[j].get("ym_i", self.shape.panels[j]["y_i"])

        for j in range(Np):
            xP, yP = control_point(j)
            t_j = tangents[j]

            # Start with freestream contribution
            V_total = U_vec.copy()

            # Add all source panel induced velocities
            for i in range(Np):
                if i == j:
                    # Self-induced tangential velocity of a source panel is zero
                    continue
                x1, y1, x2, y2 = panel_endpoints(i)
                u_vec = self._source_panel_velocity_at_point(x1, y1, x2, y2, xP, yP)
                V_total += self.sigma[i] * u_vec

            # Project onto tangent
            self.V_t[j] = np.dot(t_j, V_total)

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