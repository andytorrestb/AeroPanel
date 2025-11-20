import numpy as np

class PanelSourceSolver:
    def __init__(self, shape, U_inf, enforce_kutta=False):
        print
        self.shape = shape
        self.U_inf = U_inf
        self.enforce_kutta = enforce_kutta
        return
    
    def assemble_influence_matrix(self):
        """Assemble influence matrix using point-source approximation."""
        print("Assembling influence matrix...")
        Np = self.shape.num_panels
        self.M = np.zeros((Np, Np))

        normals = getattr(self.shape, 'N', None)
        if normals is None:
            raise ValueError("Panel normals not set on shape. Call shape.set_panel_normals() before assembling.")

        cps = self._get_control_points_array()
        S = self._get_panel_lengths_array()

        for j in range(Np):
            for i in range(Np):
                if i == j:
                    self.M[j, i] = 0.5
                    continue

                R = cps[j] - cps[i]
                r2 = np.dot(R, R)
                if r2 < 1e-18:
                    self.M[j, i] = 0.0
                else:
                    self.M[j, i] = (S[i] / (2.0 * np.pi * r2)) * np.dot(normals[j], R)

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

        if self.enforce_kutta:
            self._apply_trailing_edge_kutta_condition(U_vec)

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

        cps = self._get_control_points_array()
        S = self._get_panel_lengths_array()

        for j in range(Np):
            t_j = tangents[j]
            V_total = U_vec.copy()

            for i in range(Np):
                R = cps[j] - cps[i]
                r2 = np.dot(R, R)
                if r2 < 1e-18:
                    V_total += 0.5 * self.sigma[i] * normals[i]
                else:
                    V_total += (self.sigma[i] * S[i] / (2.0 * np.pi)) * (R / r2)

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
        # Prefer plotting vs x/c when available (suitable for airfoils)
        try:
            x = np.array([p.get("xm_i", p["x_i"]) for p in self.shape.panels])
            y = np.array([p.get("ym_i", p["y_i"]) for p in self.shape.panels])
            # Normalize by chord if attribute available
            chord = getattr(self.shape, 'chord_length', None)
            if chord is not None and chord > 0:
                x_plot = x / chord
                x_label = "x / c"
            else:
                x_plot = x
                x_label = "x"

            # Separate upper and lower surfaces based on panel ordering
            # For NACA airfoils: panels are ordered clockwise starting from TE upper surface
            # Upper surface: from TE to LE, Lower surface: from LE back toward TE
            upper_x, upper_Cp, lower_x, lower_Cp = self._separate_upper_lower_surfaces(x_plot, self.Cp, y)

            plt.figure(figsize=(10, 6))
            plt.plot(upper_x, upper_Cp, 'b-o', label='Upper Surface', linewidth=2, markersize=4)
            plt.plot(lower_x, lower_Cp, 'r-o', label='Lower Surface', linewidth=2, markersize=4)
            plt.xlabel(x_label)
            plt.ylabel("Pressure Coefficient (Cp)")
            plt.title("Pressure Coefficient Distribution")
            plt.gca().invert_yaxis()
            plt.grid(True)
            plt.legend()
            plt.tight_layout()
            plt.savefig(self.shape.output_dir + "/pressure_coefficients.png", dpi=300, bbox_inches='tight')
            plt.close()
            print("Pressure coefficient plot saved (vs x) with separated surfaces.")
        except Exception as e:
            print(f"Falling back to theta plot due to: {e}")
            theta = [panel.get("mid_angle", 0.0) for panel in self.shape.panels]
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
    
    def _separate_upper_lower_surfaces(self, x_plot, Cp, y):
        """
        Separate pressure coefficients into upper and lower surfaces for airfoils.
        
        For NACA airfoils with clockwise panel ordering:
        - Panels start at TE upper surface, go to LE along upper surface
        - Then continue from LE to TE along lower surface
        
        Returns: upper_x, upper_Cp, lower_x, lower_Cp
        """
        # Find the leading edge (minimum x position)
        min_x_idx = np.argmin(x_plot)
        
        # Split based on panel ordering and leading edge position
        # Upper surface: from start (TE) to leading edge
        upper_indices = list(range(0, min_x_idx + 1))
        # Lower surface: from leading edge to end (back toward TE)
        lower_indices = list(range(min_x_idx, len(x_plot)))
        
        # Extract upper surface data and reverse for proper plotting order (LE to TE)
        upper_x = x_plot[upper_indices][::-1]
        upper_Cp = Cp[upper_indices][::-1]
        
        # Extract lower surface data (already in LE to TE order)
        lower_x = x_plot[lower_indices]
        lower_Cp = Cp[lower_indices]
        
        return upper_x, upper_Cp, lower_x, lower_Cp
    
    def get_Cp_theta(self):
        theta = [panel["mid_angle"] for panel in self.shape.panels]
        return theta, self.Cp

    def get_Cp_x(self):
        """Return (x/c, Cp) arrays if chord_length available; otherwise (x, Cp)."""
        x = np.array([p.get("xm_i", p["x_i"]) for p in self.shape.panels])
        chord = getattr(self.shape, 'chord_length', None)
        if chord is not None and chord > 0:
            return x / chord, self.Cp
        return x, self.Cp
    
    def get_Cp_surfaces(self):
        """
        Return pressure coefficient data separated by upper and lower surfaces.
        
        Returns:
            dict: Contains 'upper' and 'lower' keys, each with 'x' and 'Cp' arrays
        """
        x = np.array([p.get("xm_i", p["x_i"]) for p in self.shape.panels])
        y = np.array([p.get("ym_i", p["y_i"]) for p in self.shape.panels])
        
        # Normalize by chord if attribute available
        chord = getattr(self.shape, 'chord_length', None)
        if chord is not None and chord > 0:
            x_plot = x / chord
        else:
            x_plot = x
            
        upper_x, upper_Cp, lower_x, lower_Cp = self._separate_upper_lower_surfaces(x_plot, self.Cp, y)
        
        return {
            'upper': {'x': upper_x, 'Cp': upper_Cp},
            'lower': {'x': lower_x, 'Cp': lower_Cp}
        }

    def _apply_trailing_edge_kutta_condition(self, U_vec):
        tangents = getattr(self.shape, 'T', None)
        if tangents is None:
            raise ValueError("Panel tangents not set. Call shape.set_panel_tangents() before enforcing Kutta condition.")

        cps = self._get_control_points_array()
        S = self._get_panel_lengths_array()
        Np = self.shape.num_panels
        if Np < 2:
            raise ValueError("Need at least two panels to enforce Kutta condition.")

        idx_upper = 0
        idx_lower = Np - 1

        row_upper = self._tangential_influence_row(idx_upper, tangents, cps, S)
        row_lower = self._tangential_influence_row(idx_lower, tangents, cps, S)

        kutta_row = row_upper - row_lower
        kutta_rhs = np.dot(U_vec, tangents[idx_lower] - tangents[idx_upper])

        self.M[-1, :] = kutta_row
        self.RHS[-1] = kutta_rhs
        print("Applied trailing-edge Kutta condition.")

    def _tangential_influence_row(self, target_idx, tangents, cps, panel_lengths):
        Np = self.shape.num_panels
        row = np.zeros(Np)
        t_vec = tangents[target_idx]
        target_cp = cps[target_idx]

        for i in range(Np):
            if i == target_idx:
                continue
            R = target_cp - cps[i]
            r2 = np.dot(R, R)
            if r2 < 1e-18:
                continue
            coef = panel_lengths[i] / (2.0 * np.pi)
            row[i] = coef * np.dot(t_vec, R) / r2

        return row

    def _get_control_points_array(self):
        cps = getattr(self.shape, 'control_points', None)
        if cps is not None and len(cps) == self.shape.num_panels:
            return np.array(cps, dtype=float)
        cp_list = []
        for panel in self.shape.panels:
            cp_list.append([
                panel.get("xm_i", panel["x_i"]),
                panel.get("ym_i", panel["y_i"])
            ])
        return np.array(cp_list, dtype=float)

    def _get_panel_lengths_array(self):
        S = getattr(self.shape, 'S', None)
        if S is None or len(S) != self.shape.num_panels:
            S = np.zeros(self.shape.num_panels)
            for i in range(self.shape.num_panels):
                j = (i + 1) % self.shape.num_panels
                dx = self.shape.panels[j]["x_i"] - self.shape.panels[i]["x_i"]
                dy = self.shape.panels[j]["y_i"] - self.shape.panels[i]["y_i"]
                S[i] = np.hypot(dx, dy)
            self.shape.S = S
        return np.array(S, dtype=float)