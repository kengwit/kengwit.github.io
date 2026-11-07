---
layout: page
title: Tensor Derivation of the Geodesic from Los Angeles to Tokyo on a Spherical Earth
description:
img: assets/img/sympy_geodesic_earth_flat.png
importance: 1
category: fun
related_publications: true
---



<iframe src="{{ '/assets/html/interactive_tensor_geodesic_globe.html' | relative_url }}" width="100%" height="500px" frameborder="0"></iframe>

# Purpose

This note derives the geodesic equations for a sphere. The derivation begins with the metric tensor of a sphere, computes the Christoffel symbols, forms the tensor geodesic equation, and then states the boundary-value problem for a path from Los Angeles to Tokyo. These results are intended for later use in the symbolic and numerical computations performed using SymPy {% cite 10.7717/peerj-cs.103 %} and SciPy.

# Coordinate System

Use 

$$
x^1=\phi,\qquad x^2=\lambda
$$ 

where $\phi$ is latitude and $\lambda$ is longitude. A sphere of radius $R$ is parametrized by

$$
\mathbf r(\phi,\lambda) =
R\begin{pmatrix}
\cos\phi\cos\lambda\\
\cos\phi\sin\lambda\\
\sin\phi
\end{pmatrix}
$$

The coordinate basis vectors are 

$$
\mathbf e_\phi = \frac{\partial\mathbf r}{\partial\phi}
=
R\begin{pmatrix}
-\sin\phi\cos\lambda\\
-\sin\phi\sin\lambda\\
\cos\phi
\end{pmatrix}
$$ 

and 

$$
\mathbf e_\lambda = \frac{\partial\mathbf r}{\partial\lambda}
=
R\begin{pmatrix}
-\cos\phi\sin\lambda\\
\cos\phi\cos\lambda\\
0
\end{pmatrix}
$$

# Metric Tensor

The metric components are 

$$
g_{ij}=\mathbf e_i\cdot\mathbf e_j
$$

Therefore, 

$$
g_{\phi\phi}=R^2,\qquad
g_{\phi\lambda}=0,\qquad
g_{\lambda\lambda}=R^2\cos^2\phi
$$ 

Thus,

$$
\boxed{
g_{ij}
=
\begin{pmatrix}
R^2&0\\
0&R^2\cos^2\phi
\end{pmatrix}
}
$$ 

and 

$$
\boxed{
ds^2=R^2d\phi^2+R^2\cos^2\phi\,d\lambda^2
}
$$

This corresponds to the SymPy construction

```
    metric = (
        R_symbol**2 * TensorProduct(dphi, dphi)
        + R_symbol**2 * cos(phi)**2
          * TensorProduct(dlam, dlam)
    )
```

The inverse metric is 

$$
\boxed{
g^{ij}
=
\begin{pmatrix}
\dfrac{1}{R^2}&0\\[6pt]
0&\dfrac{1}{R^2\cos^2\phi}
\end{pmatrix}
}
$$

# Christoffel Symbols

The Christoffel symbols of the second kind are 

$$
\boxed{
\Gamma^k_{ij}
=
\frac12g^{k\ell}
\left(
\partial_i g_{j\ell}
+\partial_j g_{i\ell}
-\partial_\ell g_{ij}
\right)
}
$$

Only one metric component depends on the coordinates:

$$
g_{\lambda\lambda}=R^2\cos^2\phi
$$ 

Therefore,

$$
\partial_\phi g_{\lambda\lambda}
=-2R^2\sin\phi\cos\phi
$$ 

and all other metric derivatives vanish.

## Calculation of $\Gamma^\phi_{\lambda\lambda}$

Set $k=\phi$ and $i=j=\lambda$. Since the metric is diagonal, only $\ell=\phi$ contributes: 

$$
\Gamma^\phi_{\lambda\lambda}
=-\frac12g^{\phi\phi}\partial_\phi g_{\lambda\lambda}
$$ 

Substitution gives 

$$
\Gamma^\phi_{\lambda\lambda}
=-\frac12\frac1{R^2}
\left(-2R^2\sin\phi\cos\phi\right)
$$ 

so 

$$
\boxed{
\Gamma^\phi_{\lambda\lambda}
=\sin\phi\cos\phi
=\frac12\sin(2\phi).
}
$$

## Calculation of $\Gamma^\lambda_{\phi\lambda}$

Only $\ell=\lambda$ contributes: 

$$
\Gamma^\lambda_{\phi\lambda}
=\frac12g^{\lambda\lambda}\partial_\phi g_{\lambda\lambda}
$$ 

Therefore,

$$
\Gamma^\lambda_{\phi\lambda}
=\frac12\frac{1}{R^2\cos^2\phi}
\left(-2R^2\sin\phi\cos\phi\right)
$$ 

which gives 

$$
\boxed{
\Gamma^\lambda_{\phi\lambda}=-\tan\phi
}
$$ 

By symmetry of the Levi--Civita connection, 

$$\boxed{
\Gamma^\lambda_{\lambda\phi}
=\Gamma^\lambda_{\phi\lambda}
=-\tan\phi
}
$$

Thus, the only nonzero Christoffel symbols are 

$$
\boxed{
\Gamma^\phi_{\lambda\lambda}=\sin\phi\cos\phi,
\qquad
\Gamma^\lambda_{\phi\lambda}
=\Gamma^\lambda_{\lambda\phi}
=-\tan\phi
}
$$

In the code, SymPy obtains them from

```
    Gamma = metric_to_Christoffel_2nd(metric)
```

# Tensor Geodesic Equation

A geodesic satisfies 

$$
\boxed{
\frac{d^2x^k}{dt^2}
+\Gamma^k_{ij}
\frac{dx^i}{dt}
\frac{dx^j}{dt}=0
}
$$ 

Equivalently, 

$$
\boxed{\nabla_{\dot\gamma}\dot\gamma=0}
$$ 

Here,

$$
\gamma(t)=
\begin{pmatrix}
\phi(t)\\
\lambda(t)
\end{pmatrix},
\qquad
\dot\gamma=\dot\phi\,\mathbf e_\phi+\dot\lambda\,\mathbf e_\lambda
$$

## Latitude Equation

For $k=\phi$, the fully expanded equation is 

$$
\ddot\phi
+\Gamma^\phi_{\phi\phi}\dot\phi^2
+\Gamma^\phi_{\phi\lambda}\dot\phi\dot\lambda
+\Gamma^\phi_{\lambda\phi}\dot\lambda\dot\phi
+\Gamma^\phi_{\lambda\lambda}\dot\lambda^2
=0
$$ 

The first three Christoffel symbols vanish, so 

$$
\ddot\phi
+\Gamma^\phi_{\lambda\lambda}\dot\lambda^2=0
$$ 

Hence, 

$$
\boxed{
\ddot\phi+\sin\phi\cos\phi\,\dot\lambda^2=0
}
$$ 

or 

$$
\boxed{
\ddot\phi=-\sin\phi\cos\phi\,\dot\lambda^2
}
$$ 

This corresponds to

```
    dy[2] = -gamma_phi_ll * lam_dot**2
```

## Longitude Equation

For $k=\lambda$, 

$$
\begin{aligned}
\ddot\lambda
&+\Gamma^\lambda_{\phi\phi}\dot\phi^2
+\Gamma^\lambda_{\phi\lambda}\dot\phi\dot\lambda\\
&+\Gamma^\lambda_{\lambda\phi}\dot\lambda\dot\phi
+\Gamma^\lambda_{\lambda\lambda}\dot\lambda^2=0
\end{aligned}
$$ 

The only nonzero terms are the two mixed terms. Thus,

$$
\ddot\lambda
-\tan\phi\,\dot\phi\dot\lambda
-\tan\phi\,\dot\lambda\dot\phi=0
$$ 

so 

$$
\boxed{
\ddot\lambda-2\tan\phi\,\dot\phi\dot\lambda=0
}
$$ 

or 

$$\boxed{
\ddot\lambda=2\tan\phi\,\dot\phi\dot\lambda
}
$$ 

This corresponds to

```
    dy[3] = -2.0 * gamma_lam_pl * phi_dot * lam_dot
```

because $\Gamma^\lambda_{\phi\lambda}=-\tan\phi$.

# First-Order System Used by the Numerical Solver

Introduce 

$$
y_1=\phi,\qquad y_2=\lambda,\qquad
y_3=\dot\phi,\qquad y_4=\dot\lambda
$$ 

Then the second-order tensor equations become 

$$
\boxed{
\begin{aligned}
\dot y_1&=y_3\\
\dot y_2&=y_4\\
\dot y_3&=-\sin y_1\cos y_1\,y_4^2\\
\dot y_4&=2\tan y_1\,y_3y_4
\end{aligned}
}
$$ 

This is the system passed to `scipy.integrate.solve_bvp`.

# Los Angeles-to-Tokyo Boundary-Value Problem

Using city-center coordinates, 

$$
\begin{aligned}
\phi_{\mathrm{LA}}&=34.0522^\circ,
&\lambda_{\mathrm{LA}}&=-118.2437^\circ,\\
\phi_{\mathrm{Tokyo}}&=35.6895^\circ,
&\lambda_{\mathrm{Tokyo}}&=139.6917^\circ
\end{aligned}
$$ 

For a continuous Pacific-going coordinate branch, Tokyo's longitude is unwrapped (subtract $360^\circ$ from Tokyo's longitude) as

$$
139.6917^\circ\equiv -220.3083^\circ
$$ 

The boundary conditions are

$$
\boxed{
\begin{aligned}
\phi(0)&=\phi_{\mathrm{LA}}\\
\lambda(0)&=\lambda_{\mathrm{LA}}\\
\phi(1)&=\phi_{\mathrm{Tokyo}}\\
\lambda(1)&=\lambda_{\mathrm{Tokyo}}
\end{aligned}
}
$$ 

The corresponding code is

```
    return np.array([
        ya[0] - phi_A,
        ya[1] - lam_A,
        yb[0] - phi_B,
        yb[1] - lam_B,
    ])
```

The initial direction is not prescribed. Instead, the boundary-value
solver finds $\dot\phi(0)$ and $\dot\lambda(0)$ such that the geodesic
beginning in Los Angeles reaches Tokyo.

# Distance Recovered from the Metric

Once $\phi(t)$ and $\lambda(t)$ are known,

$$
ds^2=R^2d\phi^2+R^2\cos^2\phi\,d\lambda^2
$$ 

Dividing by $dt^2$ gives

$$
\left(\frac{ds}{dt}\right)^2
=R^2\left(\dot\phi^2+\cos^2\phi\,\dot\lambda^2\right)
$$ 

Thus, 

$$
\boxed{
\frac{ds}{dt}
=R\sqrt{\dot\phi^2+\cos^2\phi\,\dot\lambda^2}
}
$$ 

and the total length is 

$$
\boxed{
L=R\int_0^1
\sqrt{\dot\phi^2+\cos^2\phi\,\dot\lambda^2}\,dt
}
$$ 

The implementation is

```
    metric_speed = R_KM * np.sqrt(
        phi_dot**2
        + np.cos(phi_path)**2 * lam_dot**2
    )

    distance_km = np.trapz(metric_speed, tau_plot)
```

For the numerical solution, 

$$
\boxed{L\approx8815.473\ \mathrm{km}}
$$

The path reaches a maximum latitude of approximately

$$
\boxed{47.952^\circ\ \mathrm{N}}
$$

# Calculation Chain

The complete tensor-based calculation is 

$$
\mathbf r(\phi,\lambda)
\quad\Longrightarrow\quad
g_{ij}
\quad\Longrightarrow\quad
\Gamma^k_{ij}
\quad\Longrightarrow\quad
\nabla_{\dot\gamma}\dot\gamma=0
$$ 

$$
\Longrightarrow\quad
\begin{cases}
\ddot\phi+\sin\phi\cos\phi\,\dot\lambda^2=0\\
\ddot\lambda-2\tan\phi\,\dot\phi\dot\lambda=0
\end{cases}
$$ 

$$
\Longrightarrow\quad
\text{solve the Los Angeles-to-Tokyo boundary conditions}
\quad\Longrightarrow\quad
L\approx8815.473\ \mathrm{km}
$$

This demonstrates that the route is obtained from the metric tensor and
its Levi--Civita connection rather than by directly inserting the known
great-circle distance formula $L=R\theta$.

# Results

<!-- Row 1 -->
<div class="row">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/sympy_geodesic_earth_globe.png" title="example image" class="img-fluid rounded z-depth-1" %}
    </div>
</div>
<div class="caption">
Geodesic on earth globe.
</div>

<!-- Row 2 -->
<div class="row mt-4"> <!-- Added mt-4 for a bit of spacing between rows -->
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/sympy_geodesic_earth_flat.png" title="second example image" class="img-fluid rounded z-depth-1" %}
    </div>
</div>
<div class="caption">
Geodesic on flattened earth.
</div>


<!-- 
The Python code to implement the above:

```
"""LA-to-Tokyo tensor geodesic with interactive Earth globe.

SymPy derives Christoffel symbols from the spherical metric.
SciPy solves the geodesic boundary-value problem.
Plotly writes an interactive orthographic Earth globe that can be rotated
with the mouse. A Pacific-centered flat Earth map is also written to HTML.
"""

import numpy as np
import sympy as sp
from sympy.diffgeom import Manifold, Patch, CoordSystem, TensorProduct, metric_to_Christoffel_2nd
from scipy.integrate import solve_bvp
import plotly.graph_objects as go

# ------------------------------------------------------------
# 1. Symbolic metric and Christoffel symbols
# ------------------------------------------------------------
Rsym = sp.Symbol("R", positive=True, real=True)
M = Manifold("S2", 2)
P = Patch("P", M)
S = CoordSystem("S", P, (sp.Symbol("phi", real=True), sp.Symbol("lambda", real=True)))
phi, lam = S.base_scalars()
dphi, dlam = S.base_oneforms()

metric = (
    Rsym**2 * TensorProduct(dphi, dphi)
    + Rsym**2 * sp.cos(phi)**2 * TensorProduct(dlam, dlam)
)
Gamma = metric_to_Christoffel_2nd(metric)

print("Gamma^phi_{lambda lambda} =", sp.simplify(Gamma[0, 1, 1]))
print("Gamma^lambda_{phi lambda} =", sp.simplify(Gamma[1, 0, 1]))
print("Gamma^lambda_{lambda phi} =", sp.simplify(Gamma[1, 1, 0]))

p = sp.Symbol("p", real=True)
gamma_phi_ll = sp.lambdify(p, Gamma[0, 1, 1].xreplace({phi: p}), "numpy")
gamma_lam_pl = sp.lambdify(p, Gamma[1, 0, 1].xreplace({phi: p}), "numpy")

# ------------------------------------------------------------
# 2. End points: Los Angeles and Tokyo
# ------------------------------------------------------------
R_KM = 6371.0
lat_A_deg, lon_A_deg = 34.0522, -118.2437
lat_B_deg, lon_B_deg = 35.6895, 139.6917

phi_A, lam_A = np.radians([lat_A_deg, lon_A_deg])
phi_B = np.radians(lat_B_deg)
# Unwrap Tokyo westward to select the shorter Pacific-going branch.
lam_B = np.radians(lon_B_deg - 360.0)

# ------------------------------------------------------------
# 3. Tensor geodesic equations
# ------------------------------------------------------------
def geodesic_ode(t, y):
    ph = y[0]
    ph_dot = y[2]
    lam_dot = y[3]

    return np.vstack((
        ph_dot,
        lam_dot,
        -gamma_phi_ll(ph) * lam_dot**2,
        -2.0 * gamma_lam_pl(ph) * ph_dot * lam_dot,
    ))


def boundary_conditions(ya, yb):
    return np.array([
        ya[0] - phi_A,
        ya[1] - lam_A,
        yb[0] - phi_B,
        yb[1] - lam_B,
    ])

# ------------------------------------------------------------
# 4. Solve LA -> Tokyo boundary-value problem
# ------------------------------------------------------------
t = np.linspace(0.0, 1.0, 201)
y_guess = np.zeros((4, t.size))
y_guess[0] = phi_A + (phi_B - phi_A) * t
y_guess[1] = lam_A + (lam_B - lam_A) * t
y_guess[2] = phi_B - phi_A
y_guess[3] = lam_B - lam_A

sol = solve_bvp(
    geodesic_ode,
    boundary_conditions,
    t,
    y_guess,
    tol=1.0e-10,
    max_nodes=20000,
)
if not sol.success:
    raise RuntimeError(sol.message)

t_plot = np.linspace(0.0, 1.0, 1001)
y = sol.sol(t_plot)
phi_path, lam_path, phi_dot, lam_dot = y
lat_path = np.degrees(phi_path)
lon_path = ((np.degrees(lam_path) + 180.0) % 360.0) - 180.0

metric_speed = R_KM * np.sqrt(phi_dot**2 + np.cos(phi_path)**2 * lam_dot**2)
distance_km = np.trapz(metric_speed, t_plot)
print(f"Distance from metric integral = {distance_km:.3f} km")
print(f"Maximum latitude = {lat_path.max():.3f} deg N")

# Plotly can draw an unwanted segment at the +/-180 seam if raw longitude
# jumps from -180 to +180. Insert None at that discontinuity.
route_lon, route_lat = [], []
for i, (lo, la) in enumerate(zip(lon_path, lat_path)):
    if i > 0 and abs(lo - lon_path[i - 1]) > 180.0:
        route_lon.append(None)
        route_lat.append(None)
    route_lon.append(float(lo))
    route_lat.append(float(la))

# ------------------------------------------------------------
# 5. Interactive globe with Earth geography
# ------------------------------------------------------------
fig = go.Figure()

fig.add_trace(go.Scattergeo(
    lon=route_lon,
    lat=route_lat,
    mode="lines",
    line=dict(width=4, color="crimson"),
    name="Tensor geodesic",
    hovertemplate="Lat %{lat:.2f} deg<br>Lon %{lon:.2f} deg<extra></extra>",
))

fig.add_trace(go.Scattergeo(
    lon=[lon_A_deg, lon_B_deg],
    lat=[lat_A_deg, lat_B_deg],
    mode="markers+text",
    text=["Los Angeles", "Tokyo"],
    textposition=["bottom left", "bottom right"],
    marker=dict(size=10, color=["navy", "green"], line=dict(width=1, color="white")),
    name="Endpoints",
    hovertemplate="%{text}<br>Lat %{lat:.4f} deg<br>Lon %{lon:.4f} deg<extra></extra>",
))

fig.update_geos(
    projection=dict(
        type="orthographic",
        rotation=dict(lon=-170, lat=38, roll=0),
    ),
    showland=True,
    landcolor="rgb(216,210,184)",
    showocean=True,
    oceancolor="rgb(158,202,225)",
    showcoastlines=True,
    coastlinecolor="rgb(55,55,55)",
    showcountries=True,
    countrycolor="rgb(110,110,110)",
    showlakes=True,
    lakecolor="rgb(158,202,225)",
    lataxis=dict(showgrid=True, gridcolor="rgba(100,100,100,0.25)", dtick=20),
    lonaxis=dict(showgrid=True, gridcolor="rgba(100,100,100,0.25)", dtick=30),
)

fig.update_layout(
    title=f"Los Angeles to Tokyo: tensor geodesic ({distance_km:.1f} km)",
    width=900,
    height=850,
    margin=dict(l=10, r=10, t=60, b=10),
    legend=dict(x=0.5, y=1.0, xanchor="center", yanchor="top"),
)

fig.write_html(
    "interactive_tensor_geodesic_globe.html",
    include_plotlyjs=True,
    full_html=True,
)

# ------------------------------------------------------------
# 6. Interactive flat Earth-map view for comparison
# ------------------------------------------------------------
fig_flat = go.Figure()
fig_flat.add_trace(go.Scattergeo(
    lon=route_lon,
    lat=route_lat,
    mode="lines",
    line=dict(width=4, color="crimson"),
    name="Tensor geodesic",
    hovertemplate="Lat %{lat:.2f} deg<br>Lon %{lon:.2f} deg<extra></extra>",
))
fig_flat.add_trace(go.Scattergeo(
    lon=[lon_A_deg, lon_B_deg],
    lat=[lat_A_deg, lat_B_deg],
    mode="markers+text",
    text=["Los Angeles", "Tokyo"],
    textposition=["bottom left", "bottom right"],
    marker=dict(size=10, color=["navy", "green"]),
    showlegend=False,
))
fig_flat.update_geos(
    projection=dict(type="equirectangular", rotation=dict(lon=180)),
    showland=True,
    landcolor="rgb(216,210,184)",
    showocean=True,
    oceancolor="rgb(158,202,225)",
    showcoastlines=True,
    coastlinecolor="rgb(55,55,55)",
    showcountries=True,
    countrycolor="rgb(110,110,110)",
    showlakes=True,
    lakecolor="rgb(158,202,225)",
    lataxis=dict(range=[20, 60], showgrid=True, dtick=10),
    lonaxis=dict(showgrid=True, dtick=20),
)
fig_flat.update_layout(
    title="Los Angeles to Tokyo: tensor geodesic on flat Earth map",
    width=1100,
    height=600,
    margin=dict(l=10, r=10, t=60, b=10),
)
fig_flat.write_html(
    "interactive_tensor_geodesic_flat.html",
    include_plotlyjs=True,
    full_html=True,
)

print("Wrote interactive_tensor_geodesic_globe.html")
print("Wrote interactive_tensor_geodesic_flat.html")

```
 -->
---

# Appendix: Position Vector, Coordinate Curve, and the Changing Basis

## Purpose

This appendix clarifies the distinction between the coordinate curve $\gamma(t)$ and the ambient position vector $\mathbf r(\phi,\lambda)$. It also explains why no derivative of the coordinate basis appears when one first differentiates $\gamma(t)$, even though the basis vectors vary from point to point. The basis variation appears when the tangent vector is differentiated, and intrinsically it is represented by the connection coefficients $\Gamma^k{}_{ij}$.

## Three Different Objects

The following objects should not be identified with one another.

### The coordinate tuple

A point on the coordinate patch is represented by

$$
x=(\phi,\lambda).
$$

This is a pair of coordinate values, not a vector in physical three-dimensional space.

### The curve on the manifold

A parametrized curve is a map

$$
\gamma:I\subset\mathbb R\longrightarrow S^2,
$$

which has the coordinate representation

$$
\boxed{\gamma(t)=\bigl(\phi(t),\lambda(t)\bigr).}
$$

Thus, $\gamma(t)$ identifies a point of the sphere at parameter value $t$. The displayed column $\bigl(\phi(t),\lambda(t)\bigr)^T$ contains coordinate functions. It is not the Cartesian position vector of that point.

### The embedding or position map

For a sphere of radius $R$, the embedding into $\mathbb R^3$ is

$$
\mathbf r:S^2\longrightarrow\mathbb R^3,
$$

with coordinate expression

$$
\boxed{
\mathbf r(\phi,\lambda)
=R\begin{pmatrix}
\cos\phi\cos\lambda\\
\cos\phi\sin\lambda\\
\sin\phi
\end{pmatrix}.
}
$$

The physical path in ambient space is the composition

$$
\boxed{
\mathbf c(t)=(\mathbf r\circ\gamma)(t)
=\mathbf r\bigl(\phi(t),\lambda(t)\bigr).
}
$$

Here $\mathbf c(t)\in\mathbb R^3$ is an actual Cartesian position vector. Therefore,

$$
\gamma(t)=\bigl(\phi(t),\lambda(t)\bigr)
\quad\text{and}\quad
\mathbf c(t)=\mathbf r\bigl(\phi(t),\lambda(t)\bigr)
$$

describe the same moving point, but they belong to different spaces.

## Coordinate Basis Vectors

The coordinate tangent basis is

$$
\mathbf e_\phi=\frac{\partial\mathbf r}{\partial\phi},
\qquad
\mathbf e_\lambda=\frac{\partial\mathbf r}{\partial\lambda}.
$$

Explicitly,

$$
\mathbf e_\phi
=R\begin{pmatrix}
-\sin\phi\cos\lambda\\
-\sin\phi\sin\lambda\\
\cos\phi
\end{pmatrix},
$$

and

$$
\mathbf e_\lambda
=R\begin{pmatrix}
-\cos\phi\sin\lambda\\
\cos\phi\cos\lambda\\
0
\end{pmatrix}.
$$

These basis vectors vary with $\phi$ and $\lambda$. Along the curve,

$$
\mathbf e_\phi(t)=\mathbf e_\phi\bigl(\phi(t),\lambda(t)\bigr),
\qquad
\mathbf e_\lambda(t)=\mathbf e_\lambda\bigl(\phi(t),\lambda(t)\bigr).
$$

Consequently, it would be incorrect to claim that the basis vectors are constant along a general curve. The precise statement is:

> When $\gamma(t)$ is differentiated once, the result is the tangent vector whose coordinate components are $\dot\phi$ and $\dot\lambda$. No product rule is being applied to a pre-existing basis expansion of $\gamma$.

## Why No Basis Derivative Appears in $\dot\gamma$

The expression

$$
\gamma(t)=\begin{pmatrix}\phi(t)\\\lambda(t)\end{pmatrix}
$$

is a coordinate description of a point. It is not the vector expansion

$$
\boxed{
\gamma(t)\ne\phi(t)\mathbf e_\phi+\lambda(t)\mathbf e_\lambda.
}
$$

The right-hand side would be a tangent vector, whereas $\gamma(t)$ is a point. Therefore, differentiating the coordinate tuple gives

$$
\frac{d}{dt}\bigl(\phi(t),\lambda(t)\bigr)
=\bigl(\dot\phi(t),\dot\lambda(t)\bigr).
$$

These are the components of the tangent vector in the coordinate basis:

$$
\boxed{
\dot\gamma
=\dot\phi\frac{\partial}{\partial\phi}
+\dot\lambda\frac{\partial}{\partial\lambda}.
}
$$

Equivalently, after applying the differential of the embedding $\mathbf r$, the ambient velocity is

$$
\boxed{
\dot{\mathbf c}
=\frac{d}{dt}\mathbf r\bigl(\phi(t),\lambda(t)\bigr)
=\dot\phi\mathbf e_\phi+\dot\lambda\mathbf e_\lambda.
}
$$

This is the multivariable chain rule:

$$
\frac{d\mathbf c}{dt}
=\frac{\partial\mathbf r}{\partial\phi}\frac{d\phi}{dt}
+\frac{\partial\mathbf r}{\partial\lambda}\frac{d\lambda}{dt}.
$$

At this first derivative, $\mathbf e_\phi$ and $\mathbf e_\lambda$ arise as the partial derivatives of $\mathbf r$. No derivative of those basis vectors is required yet.

## Coordinate Differentiation Does Not Freeze the Physical Basis

The notation

$$
\gamma(t)=\bigl(\phi(t),\lambda(t)\bigr)
$$

means that $\phi(t)$ and $\lambda(t)$ are the coordinate functions of the point $\gamma(t)$ in the chosen chart. Differentiating them gives the rate at which the point moves through coordinate space:

$$
\dot x^i(t)=\bigl(\dot\phi(t),\dot\lambda(t)\bigr).
$$

These numbers are the components of the tangent vector at the current point. They do not imply that the physical tangent basis is held fixed. The basis at time $t$ is

$$
\left\{
\mathbf e_\phi\bigl(\gamma(t)\bigr),
\mathbf e_\lambda\bigl(\gamma(t)\bigr)
\right\},
$$

and generally

$$
\frac{d\mathbf e_\phi}{dt}\ne0,
\qquad
\frac{d\mathbf e_\lambda}{dt}\ne0.
$$

It is useful to distinguish two operations:

1. Differentiating the coordinate map gives the component pair

   $$
   \frac{d}{dt}(\phi,\lambda)=(\dot\phi,\dot\lambda).
   $$

2. Interpreting that component pair geometrically gives the tangent vector in the basis attached to the current point:

   $$
   \dot\gamma
   =\dot\phi\mathbf e_\phi\bigl(\gamma(t)\bigr)
   +\dot\lambda\mathbf e_\lambda\bigl(\gamma(t)\bigr).
   $$

Thus, the basis is not kept constant. At each time, one uses the basis belonging to a different tangent space:

$$
\dot\gamma(t)\in T_{\gamma(t)}S^2.
$$

This also explains why ordinary subtraction of tangent vectors at different points is not intrinsically meaningful without a rule for comparing their tangent spaces. The connection supplies that rule.

## Where the Basis Derivatives Appear

Differentiate the velocity

$$
\dot{\mathbf c}=\dot\phi\mathbf e_\phi+\dot\lambda\mathbf e_\lambda.
$$

The product rule gives

$$
\begin{aligned}
\ddot{\mathbf c}
={}&\ddot\phi\mathbf e_\phi
+\ddot\lambda\mathbf e_\lambda\\
&+\dot\phi\frac{d\mathbf e_\phi}{dt}
+\dot\lambda\frac{d\mathbf e_\lambda}{dt}.
\end{aligned}
$$

The chain rule for each basis vector gives

$$
\frac{d\mathbf e_i}{dt}
=\frac{\partial\mathbf e_i}{\partial x^j}\dot x^j,
\qquad
(x^1,x^2)=(\phi,\lambda).
$$

Hence,

$$
\boxed{
\ddot{\mathbf c}
=\ddot x^k\mathbf e_k
+\dot x^i\dot x^j
\frac{\partial\mathbf e_i}{\partial x^j}.
}
$$

This is the stage at which the changing basis enters.

## How the Christoffel Symbols Account for the Moving Basis

Because the sphere is embedded in $\mathbb R^3$, the derivative of a tangent basis vector can have both tangential and normal components:

$$
\boxed{
\frac{\partial\mathbf e_i}{\partial x^j}
=\Gamma^k{}_{ij}\mathbf e_k+b_{ij}\mathbf n.
}
$$

Here,

- $\Gamma^k{}_{ij}\mathbf e_k$ is the tangential part,
- $b_{ij}\mathbf n$ is the normal part,
- $\mathbf n$ is the unit normal, and
- $b_{ij}$ is the second fundamental form.

Substitution into the ambient acceleration gives

$$
\ddot{\mathbf c}
=\left(\ddot x^k+\Gamma^k{}_{ij}\dot x^i\dot x^j\right)\mathbf e_k
+b_{ij}\dot x^i\dot x^j\mathbf n.
$$

Therefore, the tangential projection is

$$
\boxed{
(\ddot{\mathbf c})_{\mathrm{tan}}
=\left(\ddot x^k+\Gamma^k{}_{ij}\dot x^i\dot x^j\right)\mathbf e_k.
}
$$

The covariant acceleration is this tangential part:

$$
\boxed{
\nabla_{\dot\gamma}\dot\gamma
=\left(\ddot x^k+\Gamma^k{}_{ij}\dot x^i\dot x^j\right)
\frac{\partial}{\partial x^k}.
}
$$

This formula can be read as

$$
\boxed{
\text{covariant change}
=
\text{component change}
+
\text{change caused by the moving basis}.
}
$$

The term $\ddot x^k$ differentiates the coordinate components as ordinary scalar functions. The Christoffel term restores the effect of the changing basis. Thus, the Christoffel symbols do not assert that the basis is constant; they are precisely the quantities that account for the fact that it is not constant.

## Geodesic Equation

A geodesic has zero tangential acceleration:

$$
\boxed{\nabla_{\dot\gamma}\dot\gamma=0.}
$$

In coordinates,

$$
\boxed{
\ddot x^k+\Gamma^k{}_{ij}\dot x^i\dot x^j=0.
}
$$

This does not generally imply $\ddot{\mathbf c}=0$. Instead, it implies that $\ddot{\mathbf c}$ has no tangential component. For a geodesic on a sphere, the ambient acceleration is purely normal to the surface.

For latitude-longitude coordinates, the nonzero Christoffel symbols are

$$
\Gamma^\phi{}_{\lambda\lambda}=\sin\phi\cos\phi,
\qquad
\Gamma^\lambda{}_{\phi\lambda}
=\Gamma^\lambda{}_{\lambda\phi}=-\tan\phi.
$$

Therefore,

$$
\boxed{
\ddot\phi+\sin\phi\cos\phi\dot\lambda^2=0,
}
$$

and

$$
\boxed{
\ddot\lambda-2\tan\phi\dot\phi\dot\lambda=0.
}
$$

The additional terms are not artificial corrections. They are precisely the terms required because the coordinate basis changes along the path.

## Consistency Check in Cartesian Coordinates

In flat Euclidean space with a global Cartesian basis,

$$
\frac{d\mathbf e_i}{dt}=0
\quad\Longrightarrow\quad
\Gamma^k{}_{ij}=0.
$$

The geodesic equation then reduces to

$$
\ddot x^k=0,
$$

which is the ordinary straight-line equation. In spherical coordinates, the basis changes from point to point, so the Christoffel symbols are generally nonzero. The resulting extra terms are the coordinate expression of going straight while using a moving, non-Cartesian basis.

## A General Vector Field Along the Curve

For a general tangent vector field

$$
V(t)=V^i(t)\mathbf e_i\bigl(\gamma(t)\bigr),
$$

its ordinary ambient derivative is

$$
\frac{dV}{dt}
=\dot V^k\mathbf e_k
+V^i\dot x^j\frac{\partial\mathbf e_i}{\partial x^j}.
$$

Taking the tangential part gives

$$
\boxed{
\nabla_{\dot\gamma}V
=\left(\dot V^k+\Gamma^k{}_{ij}\dot x^iV^j\right)\mathbf e_k.
}
$$

For parallel transport,

$$
\boxed{\nabla_{\dot\gamma}V=0.}
$$

Thus, a vector can have changing coordinate or Cartesian components while having zero intrinsic tangential change.

## Summary of the Distinction

| Object | Meaning |
|---|---|
| $\gamma(t)$ | A point-valued curve on the manifold, represented by the coordinate pair $(\phi(t),\lambda(t))$. |
| $\mathbf r(\phi,\lambda)$ | The embedding or position map from spherical coordinates into $\mathbb R^3$. |
| $\mathbf c(t)=\mathbf r\circ\gamma$ | The Cartesian position vector of the moving point. |
| $\dot\gamma$ | The abstract tangent vector with coordinate components $(\dot\phi,\dot\lambda)$. |
| $\dot{\mathbf c}$ | The ambient representation $\dot\phi\mathbf e_\phi+\dot\lambda\mathbf e_\lambda$. |
| $\Gamma^k{}_{ij}$ | The tangential part of the changing coordinate basis, used when differentiating tangent vectors. |

The key point is

$$
\boxed{
\gamma(t)\ne\phi(t)\mathbf e_\phi+\lambda(t)\mathbf e_\lambda.
}
$$

Therefore, differentiating $\gamma(t)$ does not call for a product rule on basis vectors. The basis variation first appears when differentiating the tangent vector $\dot\gamma$, and the covariant derivative records that variation through the Christoffel symbols.
