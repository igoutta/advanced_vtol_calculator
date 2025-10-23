import marimo

__generated_with = "0.17.0"
app = marimo.App(
    width="medium",
    auto_download=["ipynb"],
    sql_output="lazy-polars",
)


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _():
    from math import cos, degrees, radians, sin, sqrt
    return cos, degrees, radians, sin, sqrt


@app.cell
def _():
    import numpy as np
    return (np,)


@app.cell
def _():
    import sympy as sp
    from sympy import Eq, linsolve, solve, symbols
    return


@app.cell
def _(mo):
    mo.md(r"""## First iteration""")
    return


@app.cell
def _():
    # Gravitational force
    m = 25.5  # kg o kgf
    g = 9.81  # m/s^2
    G = m * g
    return G, g


@app.cell
def _():
    # Pitch
    th = 10
    # Roll
    ph = 0
    return ph, th


@app.cell
def _(G, g, ph, radians, sin, sqrt, th):
    # First iteration
    if abs(th) > 90 or abs(ph) > 90:
        raise AttributeError("Any angle can be over 90 degrees")

    firstT = G / sqrt(1 - sin(radians(th)) ** 2 - sin(radians(ph)) ** 2)
    print(f"The thrust per motor should be {firstT / (g * 4):.3f} kgf")
    return (firstT,)


@app.cell
def _(firstT, g, ph, radians, sin, th):
    # Max drag forces it could take
    f1Fx = firstT * sin(radians(th))
    f1Fy = firstT * sin(radians(ph))

    print(f"In X axis, the UAV exceed {f1Fx:.3f} N or {f1Fx / g:.3f} kgf")
    print(f"In Y axis, the UAV exceed {f1Fy:.3f} N or {f1Fy / g:.3f} kgf")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    # Equations

    $$
    i_x = T_1\sin{\theta} - T_2\sin{\theta} + T_3\sin{\theta} - T_4\sin{\theta}
    $$

    $$
    i_y = - T_1\sin{\phi} + T_2\sin{\phi} + T_3\sin{\phi} - T_4\sin{\phi}
    $$

    $$
    i_z + G = (T_1 + T_2 + T_3 + T_4)\sqrt{1 - \sin^2{\theta} - \sin^2{\phi}}
    $$


    $$
    \begin{align}
    \vec{M_i} &= c\vec{T_1} + \vec{r_1}\times\vec{T_1} + \\
              &= c\vec{T_2} + \vec{r_2}\times\vec{T_2} + \\
              &= c\vec{T_3} + \vec{r_3}\times\vec{T_3} + \\
              &= c\vec{T_4} + \vec{r_4}\times\vec{T_4}
    \end{align}
    $$
    """
    )
    return


@app.cell
def _(mo):
    mo.md(
        r"""
    Rotation Matrix around X(roll) and Y(pitch) for radius r : $\vec{R} = \Phi\times\Theta\times\vec{r}$

    $$
    \begin{pmatrix}
    1 & 0 & 0\\
    0 & \cos(\phi) & -\sin(\phi) \\
    0 & \sin(\phi) & \cos(\phi)
    \end{pmatrix}
    \begin{pmatrix}
    \cos(\theta) & 0 & \sin(\theta) \\
    0 & 1 & 0\\
    -\sin(\theta) & 0 & \cos(\theta))
    \end{pmatrix}
    \begin{pmatrix}
    r\cos(\alpha)\\
    r\sin(\alpha)\\
    0
    \end{pmatrix}
    =
    \begin{pmatrix}
    r\cos(\alpha)\cos(\theta) \\
    r[\sin(\alpha)\cos(\phi) + \cos(\alpha)\sin(\theta)\sin(\phi)] \\
    r[\sin(\alpha)\sin(\phi) - \cos(\alpha)\sin(\theta)\cos(\phi)]
    \end{pmatrix}
    $$
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Inputs""")
    return


@app.cell
def _():
    i_x = 5  # kgf
    i_y = 4  # kgf
    i_z = 0  # kgf
    m_ix = 0  # kgf-m
    m_iy = 0  # kgf-m
    m_iz = 0  # kgf-m
    return i_x, i_y, i_z


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Outputs""")
    return


@app.cell
def _(G, i_x, i_y, i_z):
    # EF(x)
    T_x = i_x * 9.80665
    # EF(y)
    T_y = i_y * 9.80665
    # EF(z)
    T_z = G + i_z
    return T_x, T_y, T_z


@app.cell
def _(T_x, T_y, T_z, degrees, np, sqrt):
    norm_T = sqrt(T_x**2 + T_y**2 + T_z**2)

    theta = degrees(np.asin(T_x / norm_T))
    phi = degrees(np.asin(T_y / norm_T))
    print(f"Pitch: {theta}° | Roll: {phi}°")
    return phi, theta


@app.cell
def _(np):
    T = dict()
    T[1] = np.array([[1, 1, 1]]).T
    T[2] = np.ones((3, 1))
    T[3] = np.ones((3, 1))
    T[4] = np.ones((3, 1))
    T
    return (T,)


@app.cell
def _(cos, np, radians, sin):
    def get_radius_matrix(r, alpha, theta, phi):
        return np.array(
            [
                [
                    r * cos(radians(theta)) * cos(radians(alpha)),
                    r
                    * (
                        sin(radians(alpha)) * cos(radians(phi))
                        + cos(radians(alpha)) * sin(radians(theta)) * sin(radians(phi))
                    ),
                    r
                    * (
                        sin(radians(alpha)) * sin(radians(phi))
                        + cos(radians(alpha)) * sin(radians(theta)) * cos(radians(phi))
                    ),
                ]
            ]
        )
    return (get_radius_matrix,)


@app.cell
def _(get_radius_matrix, phi, theta):
    R = dict()
    R[1] = get_radius_matrix(0.8, 45, theta, phi)
    R[2] = get_radius_matrix(0.8, 45 + 180, theta, phi)
    R[3] = get_radius_matrix(0.8, 45 + 90, theta, phi)
    R[4] = get_radius_matrix(0.8, 360 - 45, theta, phi)
    R
    return (R,)


@app.cell
def _():
    cM = dict()
    cM[1] = 0
    cM[2] = 0
    cM[3] = 0
    cM[4] = 0
    return (cM,)


@app.cell
def _(R, T, cM, np):
    # EM
    M = dict()
    M[1] = cM[1] * np.linalg.norm(T[1]) + T[1] @ R[1]
    M[2] = cM[2] * np.linalg.norm(T[2]) + T[2] @ R[2]
    M[3] = cM[3] * np.linalg.norm(T[3]) + T[3] @ R[3]
    M[4] = cM[4] * np.linalg.norm(T[4]) + T[4] @ R[4]
    M
    return (M,)


@app.cell
def _(M, np):
    sumatoria_momentos = np.sum(list(M.values())).round()
    print(sumatoria_momentos)
    return


if __name__ == "__main__":
    app.run()
