import math
from typing import Tuple, List

import numpy as np
import matplotlib.pyplot as plt

from geometry.Shape import Shape


class NACA(Shape):
    """NACA 4-digit airfoil geometry and panelization.

    This class generates a closed polygon of panels around a 4-digit NACA airfoil
    using either uniform or cosine spacing. The panel ordering is clockwise,
    starting at the trailing edge along the upper surface to the leading edge,
    then along the lower surface back toward (but not including) the trailing edge.

    Notes:
    - Panels are stored similar to Circle: a list of dicts with keys x_i, y_i and
      also xm_i, ym_i (panel midpoints) once set_panels() is called.
    - Tangents and normals are computed per panel; for clockwise ordering, the
      outward normal is the +90° (CCW) rotation of the unit tangent.
    - We use the 0.1015 coefficient in the thickness polynomial so the trailing
      edge has a small finite thickness, avoiding a zero-length panel at the seam.
    """

    def __init__(self, designation: str, chord_length: float = 1.0):
        self.designation = designation.strip()
        self.chord_length = chord_length

        # Populated by set_panels
        self.num_panels: int = 0
        self.panels: List[dict] = []

    # ----------------------------
    # Public API
    # ----------------------------
    def set_panels(
        self,
        num_panels: int,
        spacing: str = "cosine",
        angle_of_attack_deg: float = 0.0,
        origin: Tuple[float, float] = (0.0, 0.0),
    ) -> None:
        """Discretize the NACA 4-digit airfoil into straight-line panels.

        Parameters:
        - num_panels: total number of panels (must be an even integer >= 4).
        - spacing: 'cosine' (default) or 'uniform' chordwise spacing per surface.
        - angle_of_attack_deg: rotate the airfoil by AoA (degrees) about origin.
        - origin: rotation/translation origin for the airfoil (x0, y0).
        """
        if num_panels < 4 or (num_panels % 2) != 0:
            raise ValueError("num_panels must be an even integer >= 4")

        m, p, t = self._parse_naca_4(self.designation)

        # Number of x-stations per surface including both LE and TE
        K = num_panels // 2 + 1
        x01 = self._x_distribution(K, spacing)

        # Base chord coordinates (before rotation/translation)
        x = x01 * self.chord_length

        # Thickness distribution (with finite TE thickness via 0.1015 coefficient)
        yt = self._thickness_dist(x01, t) * self.chord_length

        # Camber line and slope
        yc, dyc_dx = self._camber_and_slope(x01, m, p)

        # Surface coordinates about camber line
        theta = np.arctan(dyc_dx)
        xu = x - yt * np.sin(theta)
        yu = yc * self.chord_length + yt * np.cos(theta)
        xl = x + yt * np.sin(theta)
        yl = yc * self.chord_length - yt * np.cos(theta)

        # Build clockwise boundary point list:
        # - Start at TE (upper), go to LE along upper: reverse upper
        # - Then from LE toward TE along lower, excluding LE and TE
        upper_pts = np.stack((xu, yu), axis=1)[::-1]  # TE->LE
        lower_pts = np.stack((xl, yl), axis=1)        # LE->TE

        # Exclude LE and TE on lower to avoid duplicates; we'll close the loop via wrap-around
        lower_slice = lower_pts[1:-1] if lower_pts.shape[0] > 2 else np.empty((0, 2))

        pts = np.vstack((upper_pts, lower_slice))
        if pts.shape[0] != num_panels:
            # In case spacing/rounding caused mismatch, trim or pad conservatively
            if pts.shape[0] > num_panels:
                pts = pts[:num_panels]
            else:
                # Pad by repeating last point slightly jittered in x to keep topology
                deficit = num_panels - pts.shape[0]
                if pts.shape[0] == 0:
                    raise RuntimeError("Failed to generate any panel points for NACA airfoil.")
                last = pts[-1]
                jitter = 1e-9
                pads = np.array([last + [jitter * (i + 1), 0.0] for i in range(deficit)])
                pts = np.vstack((pts, pads))

        # Apply angle of attack rotation and translation about origin
        aoa = math.radians(angle_of_attack_deg)
        cos_a, sin_a = math.cos(aoa), math.sin(aoa)
        x0, y0 = origin
        R = np.array([[cos_a, -sin_a], [sin_a, cos_a]])
        pts_rot = (R @ pts.T).T + np.array([[x0, y0]])

        # Populate panel list with boundary start points and midpoints
        self.num_panels = num_panels
        self.panels = []
        for i in range(num_panels):
            j = (i + 1) % num_panels
            x_i, y_i = pts_rot[i]
            x_j, y_j = pts_rot[j]
            xm, ym = 0.5 * (x_i + x_j), 0.5 * (y_i + y_j)

            # Angle of panel orientation (for optional plotting compatibility)
            beta = math.degrees(math.atan2(y_j - y_i, x_j - x_i))

            self.panels.append(
                {
                    "x_i": float(x_i),
                    "y_i": float(y_i),
                    "xm_i": float(xm),
                    "ym_i": float(ym),
                    "mid_angle": beta,
                }
            )

        # Compute dependent geometric arrays
        self._set_panel_vectors_from_points()

    def set_panels_cosine(
        self,
        num_panels: int,
        angle_of_attack_deg: float = 0.0,
        origin: Tuple[float, float] = (0.0, 0.0),
    ) -> None:
        """Convenience wrapper for set_panels(..., spacing='cosine')."""
        self.set_panels(
            num_panels=num_panels,
            spacing="cosine",
            angle_of_attack_deg=angle_of_attack_deg,
            origin=origin,
        )

    def set_panels_uniform(
        self,
        num_panels: int,
        angle_of_attack_deg: float = 0.0,
        origin: Tuple[float, float] = (0.0, 0.0),
    ) -> None:
        """Convenience wrapper for set_panels(..., spacing='uniform')."""
        self.set_panels(
            num_panels=num_panels,
            spacing="uniform",
            angle_of_attack_deg=angle_of_attack_deg,
            origin=origin,
        )

    def set_panel_normals(self) -> None:
        """Compute outward unit normals for each panel from boundary points.

        For clockwise panel ordering, n_out = rotate_CCW_90(t_hat).
        """
        if not self.panels or self.num_panels == 0:
            raise ValueError("Panels not set. Call set_panels() first.")
        # Determine boundary orientation by signed area
        area = self._signed_area_from_panels()
        # If area < 0 => CW; area > 0 => CCW
        rotate = "ccw" if area < 0 else "cw"

        N = np.zeros((self.num_panels, 2))
        for i in range(self.num_panels):
            j = (i + 1) % self.num_panels
            xi, yi = self.panels[i]["x_i"], self.panels[i]["y_i"]
            xj, yj = self.panels[j]["x_i"], self.panels[j]["y_i"]
            dx, dy = xj - xi, yj - yi
            L = math.hypot(dx, dy)
            if L == 0.0:
                N[i] = np.array([0.0, 0.0])
                continue
            tx, ty = dx / L, dy / L
            if rotate == "ccw":
                # Outward = +90° rotation when ordering is CW
                nx, ny = -ty, tx
            else:
                # Outward = -90° rotation when ordering is CCW
                nx, ny = ty, -tx
            N[i] = np.array([nx, ny])
        self.N = N

    def set_panel_tangents(self) -> None:
        """Compute unit tangent vectors for each panel from boundary points."""
        if not self.panels or self.num_panels == 0:
            raise ValueError("Panels not set. Call set_panels() first.")

        T = np.zeros((self.num_panels, 2))
        for i in range(self.num_panels):
            j = (i + 1) % self.num_panels
            xi, yi = self.panels[i]["x_i"], self.panels[i]["y_i"]
            xj, yj = self.panels[j]["x_i"], self.panels[j]["y_i"]
            dx, dy = xj - xi, yj - yi
            L = math.hypot(dx, dy)
            if L == 0.0:
                T[i] = np.array([0.0, 0.0])
            else:
                T[i] = np.array([dx / L, dy / L])
        self.T = T

    def print_info(self):
        m, p, t = self._parse_naca_4(self.designation)
        print(
            f"NACA {self.designation}: m={m:.3f}, p={p:.3f}, t={t:.3f}, chord={self.chord_length}"
        )
        print(f"Number of panels: {self.num_panels}")
        if not self.panels:
            print("Panels not initialized. Call set_panels(num_panels) first.")
            return
        for i, panel in enumerate(self.panels):
            print(
                f"Panel {i}: x_i={panel['x_i']:.5f}, y_i={panel['y_i']:.5f}, "
                f"xm_i={panel['xm_i']:.5f}, ym_i={panel['ym_i']:.5f}, mid_angle={panel['mid_angle']:.2f}°"
            )

    # ----------------------------
    # Plotting (axis-scaled to airfoil extents)
    # ----------------------------
    def _axes_limits(self, margin_ratio: float = 0.1):
        xs = np.array([p["x_i"] for p in self.panels])
        ys = np.array([p["y_i"] for p in self.panels])
        xmin, xmax = xs.min(), xs.max()
        ymin, ymax = ys.min(), ys.max()
        xr = xmax - xmin
        yr = ymax - ymin
        mrgx = xr * margin_ratio if xr > 0 else 0.1
        mrgy = yr * margin_ratio if yr > 0 else 0.1
        return (xmin - mrgx, xmax + mrgx, ymin - mrgy, ymax + mrgy)

    def plot_airfoil(self):
        if not self.panels:
            raise ValueError("Panels not set. Call set_panels() first.")

        fig, ax = plt.subplots()

        x_i = np.array([p["x_i"] for p in self.panels])
        y_i = np.array([p["y_i"] for p in self.panels])
        # Close loop for line plot
        x_plot = np.r_[x_i, x_i[0]]
        y_plot = np.r_[y_i, y_i[0]]

        ax.plot(x_plot, y_plot, 'r-', label='Panel Edges')
        ax.plot(x_i, y_i, 'ro', label='Panel Start Points')

        ax.set_aspect('equal', 'box')
        xmin, xmax, ymin, ymax = self._axes_limits()
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_title(f'NACA {self.designation} with Panels')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.grid(True)
        ax.legend(loc='center', bbox_to_anchor=(0.5, -1.85), frameon=True)
        plt.savefig(self.output_dir + "/airfoil_with_panels.png")
        plt.close(fig)

    def plot_normals(self, scale: float = 1.0):
        if not hasattr(self, 'N'):
            raise ValueError("Panel normals not set. Call set_panel_normals() first.")

        fig, ax = plt.subplots()
        ax.set_title("Panel Normals")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")

        for i in range(self.num_panels):
            x_i = self.panels[i]["xm_i"]
            y_i = self.panels[i]["ym_i"]
            nx, ny = self.N[i]
            ax.arrow(x_i, y_i, scale * nx, scale * ny,
                     head_width=0.02 * self.chord_length,
                     head_length=0.04 * self.chord_length,
                     fc='b', ec='b')

        # Draw outline too for reference
        xs = [p["x_i"] for p in self.panels] + [self.panels[0]["x_i"]]
        ys = [p["y_i"] for p in self.panels] + [self.panels[0]["y_i"]]
        ax.plot(xs, ys, 'k-', linewidth=1.0, alpha=0.5)

        ax.set_aspect('equal', 'box')
        xmin, xmax, ymin, ymax = self._axes_limits()
        ax.set_xlim(xmin - 0.25, xmax + 0.25)
        ax.set_ylim(ymin - 0.5, ymax + 0.5)
        ax.grid(True)
        plt.savefig(self.output_dir + "/panel_normals.png")
        plt.close(fig)

    def plot_tangents(self, scale: float = 1.0):
        if not hasattr(self, 'T'):
            raise ValueError("Panel tangents not set. Call set_panel_tangents() first.")

        fig, ax = plt.subplots()
        ax.set_title("Panel Tangents")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")

        for i in range(self.num_panels):
            x_i = self.panels[i]["xm_i"]
            y_i = self.panels[i]["ym_i"]
            tx, ty = self.T[i]
            ax.arrow(x_i, y_i, scale * tx, scale * ty,
                     head_width=0.02 * self.chord_length,
                     head_length=0.04 * self.chord_length,
                     fc='g', ec='g')

        xs = [p["x_i"] for p in self.panels] + [self.panels[0]["x_i"]]
        ys = [p["y_i"] for p in self.panels] + [self.panels[0]["y_i"]]
        ax.plot(xs, ys, 'k-', linewidth=1.0, alpha=0.5)

        ax.set_aspect('equal', 'box')
        xmin, xmax, ymin, ymax = self._axes_limits()
        ax.set_xlim(xmin - 0.25, xmax + 0.25)
        ax.set_ylim(ymin - 0.5, ymax + 0.5)
        ax.grid(True)
        plt.savefig(self.output_dir + "/panel_tangents.png")
        plt.close(fig)

    def plot_control_points(self):
        if not hasattr(self, 'control_points'):
            raise ValueError("Control points not set. Call set_panels() first (or set_control_points()).")

        fig, ax = plt.subplots()
        ax.set_title("Panel Control Points")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")

        ax.plot(self.control_points[:, 0], self.control_points[:, 1], 'ms', label='Control Points')

        xs = [p["x_i"] for p in self.panels] + [self.panels[0]["x_i"]]
        ys = [p["y_i"] for p in self.panels] + [self.panels[0]["y_i"]]
        ax.plot(xs, ys, 'k-', linewidth=1.0, alpha=0.5, label='Airfoil Outline')

        ax.set_aspect('equal', 'box')
        xmin, xmax, ymin, ymax = self._axes_limits()
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.grid(True)
        ax.legend(loc='center', bbox_to_anchor=(0.5, -1.85), frameon=True)
        plt.savefig(self.output_dir + "/panel_control_points.png")
        plt.close(fig)

    # ----------------------------
    # Internals
    # ----------------------------
    @staticmethod
    def _parse_naca_4(designation: str) -> Tuple[float, float, float]:
        """Parse a 4-digit NACA designation into (m, p, t).

        Accepts forms like '2412' or 'NACA2412' or 'NACA 2412'.
        Returns:
          m = maximum camber (fraction of chord)
          p = location of max camber (fraction of chord)
          t = maximum thickness (fraction of chord)
        """
        s = designation.upper().replace("NACA", "").replace(" ", "").strip()
        if len(s) != 4 or not s.isdigit():
            raise ValueError(f"Invalid 4-digit NACA designation: '{designation}'")
        m = int(s[0]) / 100.0
        p = int(s[1]) / 10.0
        t = int(s[2:]) / 100.0
        return m, p, t

    @staticmethod
    def _x_distribution(K: int, spacing: str) -> np.ndarray:
        """Return K points x/c in [0, 1], including both 0 and 1."""
        if spacing.lower() == "cosine":
            i = np.arange(K)
            return 0.5 * (1.0 - np.cos(np.pi * i / (K - 1)))
        elif spacing.lower() == "uniform":
            return np.linspace(0.0, 1.0, K)
        else:
            raise ValueError("spacing must be 'cosine' or 'uniform'")

    @staticmethod
    def _thickness_dist(x: np.ndarray, t: float) -> np.ndarray:
        """NACA 4-digit thickness distribution y_t/c for x in [0,1].

        Uses 0.1015 for the x^4 coefficient (finite TE thickness).
        """
        return 5.0 * t * (
            0.2969 * np.sqrt(np.clip(x, 0.0, 1.0))
            - 0.1260 * x
            - 0.3516 * x**2
            + 0.2843 * x**3
            - 0.1015 * x**4
        )

    @staticmethod
    def _camber_and_slope(x: np.ndarray, m: float, p: float) -> Tuple[np.ndarray, np.ndarray]:
        """Return camber y_c/c and slope dyc/dx for x in [0,1]."""
        if m == 0.0 or p == 0.0:
            yc = np.zeros_like(x)
            dyc = np.zeros_like(x)
            return yc, dyc

        yc = np.zeros_like(x)
        dyc = np.zeros_like(x)
        for idx, xi in enumerate(x):
            if xi < p:
                yc[idx] = (m / (p**2)) * (2 * p * xi - xi**2)
                dyc[idx] = (2 * m / (p**2)) * (p - xi)
            else:
                yc[idx] = (m / ((1 - p) ** 2)) * ((1 - 2 * p) + 2 * p * xi - xi**2)
                dyc[idx] = (2 * m / ((1 - p) ** 2)) * (p - xi)
        return yc, dyc

    def _set_panel_vectors_from_points(self) -> None:
        """Set panel lengths (S), normals (N), tangents (T), and control points."""
        # Panel lengths
        S = np.zeros(self.num_panels)
        N = np.zeros((self.num_panels, 2))
        T = np.zeros((self.num_panels, 2))
        cps = []

        # Determine boundary orientation by signed area
        area = self._signed_area_from_panels()
        rotate = "ccw" if area < 0 else "cw"
        for i in range(self.num_panels):
            j = (i + 1) % self.num_panels
            xi, yi = self.panels[i]["x_i"], self.panels[i]["y_i"]
            xj, yj = self.panels[j]["x_i"], self.panels[j]["y_i"]
            dx, dy = xj - xi, yj - yi
            L = math.hypot(dx, dy)
            S[i] = L
            if L == 0.0:
                T[i] = np.array([0.0, 0.0])
                N[i] = np.array([0.0, 0.0])
            else:
                tx, ty = dx / L, dy / L
                T[i] = np.array([tx, ty])
                if rotate == "ccw":
                    N[i] = np.array([-ty, tx])  # outward for CW ordering
                else:
                    N[i] = np.array([ty, -tx])  # outward for CCW ordering
            # 3/4 control point like Shape.set_control_points
            cx = xi + 0.75 * dx
            cy = yi + 0.75 * dy
            cps.append((cx, cy))

        self.S = S
        self.N = N
        self.T = T
        self.control_points = np.array(cps)

    def _signed_area_from_panels(self) -> float:
        """Compute polygon signed area from panel start points (shoelace formula).

        > 0 => CCW, < 0 => CW
        """
        x = np.array([p["x_i"] for p in self.panels])
        y = np.array([p["y_i"] for p in self.panels])
        x_next = np.roll(x, -1)
        y_next = np.roll(y, -1)
        return 0.5 * np.sum(x * y_next - x_next * y)