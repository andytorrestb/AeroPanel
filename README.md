# AeroPanel - Pseudocode

This pseudocode outlines the key steps to implement a **2D external flow solver** using the **constant-strength source panel method** (no Kutta condition, no vortices).

---

## Inputs
- Body geometry (e.g., circle, NACA 0012)
- Number of panels `N`
- Freestream velocity `U∞` and angle of attack `α`\

---

## Step-by-Step Algorithm

### 1. Generate Panel Geometry
- Discretize surface into `N` panels with endpoints `(x_i, y_i)` (clockwise order)
- For each panel:
  - Compute panel length `Δs_i`
  - Compute panel unit tangent `t̂_i`
  - Compute outward unit normal `n̂_i`
  - Place control point at 3/4 of the panel

---

### 2. Compute Influence Coefficients
- For each control point `j`:
  - For each panel `i`:
    - Calculate the source panel influence `u^(S)(x_j; σ_i = 1)`
    - Project onto the normal:  
      `A_ji^(S) = n̂_j ⋅ u^(S)(x_j; σ_i = 1)`
- Construct the matrix `A` of size `N x N`

---

### 3. Assemble Right-Hand Side
- For each control point `j`:
  - Compute `RHS_j = - n̂_j ⋅ U∞`

---

### 4. Solve Linear System
- Solve `A * σ = RHS` for `σ_i` (panel source strengths)

---

### 5. Compute Tangential Velocity
- For each control point `j`:
  - Compute tangential velocity:
    `u_t,j = t̂_j ⋅ ( U∞ + ∑ u^(S)(x_j; σ_i) )`

---

### 6. Compute Pressure Coefficient
- For each panel `j`:
  - `Cp_j = 1 - (u_t,j / U∞)^2`

---

### 7. Integrate for Forces and Moments
- Use numerical integration (e.g., trapezoidal rule) to compute:
  - `Cl` – Lift coefficient
  - `Cm` – Moment coefficient (about quarter chord)

---

### 8. Plot and Output
- Plot:
  - Geometry with normals and control points
  - `Cp` distribution
- Output:
  - `σ_i`, `u_t`, `Cp`, `Cl`, `Cm`

---

## Test Cases
- **Case 1: Circular Cylinder**  
  Compare `Cp(θ)` with analytic solution:  
  `Cp(θ) = 1 - 4 sin²(θ)`

- **Case 2: NACA 0012**  
  Run for `α = 0°` and `5°`, analyze nonphysical trailing edge behavior due to lack of Kutta condition.

- **Case 3: XFOIL Comparison**  
  Compare results (Cp, Cl, Cm) with XFOIL output at same angle of attack.

---

## Notes
- Use **cosine spacing** for better clustering near leading/trailing edges
- Ensure **self-influence** is correctly computed
- Use **double precision** to handle matrix conditioning as `N` increases

---

## Data-Flow Diagram

```mermaid
flowchart TD
  U["Freestream (U∞, α)"]
  G["Geometry (shape, N)"]
  Geo["Panels + control pts\nΔs, t̂, n̂"]
  A["Influence matrix A\n(normal velocity)"]
  RHS["RHS = -n̂·U∞"]
  Solve["Solve for σ"]
  Ut["Tangential velocity u_t\n(t̂·(U∞ + induced))"]
  Cp["Pressure Cp"]
  Forces["Integrate → Cl, Cm"]
  Output["Plots + tables"]

  U --> RHS
  G --> Geo
  Geo --> A
  Geo --> RHS
  Geo --> Ut
  A --> Solve
  RHS --> Solve
  Solve --> Ut
  Ut --> Cp
  Cp --> Forces
  Cp --> Output
  Forces --> Output

  subgraph Legend
    L1["Δs: panel length"]
    L2["t̂ / n̂: tangent / normal unit vectors"]
    L3["σ: source strength"]
    L4["u_t: tangential velocity"]
    L5["Cp: pressure coefficient"]
    L6["Cl, Cm: lift / moment coefficients"]
  end
```
