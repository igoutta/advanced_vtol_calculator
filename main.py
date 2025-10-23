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
    from math import degrees, radians
    return (radians,)


@app.cell
def _():
    import numpy as np
    return (np,)


@app.cell
def _():
    import sympy as sp
    from sympy import Eq, linsolve, solve, symbols
    from sympy import cos, sin, sqrt
    return Eq, cos, sin, solve, sp, sqrt, symbols


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
def _(G, g, np, ph, radians, th):
    # First iteration
    if abs(th) > 90 or abs(ph) > 90:
        raise AttributeError("Any angle can be over 90 degrees")

    firstT = G / np.sqrt(1 - np.sin(radians(th)) ** 2 - np.sin(radians(ph)) ** 2)
    print(f"The thrust per motor should be {firstT / (g * 4):.3f} kgf")
    return (firstT,)


@app.cell
def _(firstT, g, np, ph, radians, th):
    # Max drag forces it could take
    f1Fx = firstT * np.sin(radians(th))
    f1Fy = firstT * np.sin(radians(ph))

    print(f"In X axis, the UAV exceed {f1Fx:.3f} N or {f1Fx / g:.3f} kgf")
    print(f"In Y axis, the UAV exceed {f1Fy:.3f} N or {f1Fy / g:.3f} kgf")
    return


@app.cell
def _(mo):
    mo.md(
        r"""
    ## Physics behind

    $$
    \begin{equation}\vec{M} = \vec{r}\times\vec{F}\end{equation}
    $$
    """
    )
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
    $$
    T_n = 
    \begin{pmatrix}
    ||T_n||\sin{\theta} \\
    ||T_n||\sin{\phi} \\
    ||T_n||\sqrt{1 - \sin^2{\theta} - \sin^2{\phi}}
    \end{pmatrix}
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
    0 & \cos{\phi} & -\sin{\phi} \\
    0 & \sin{\phi} & \cos{\phi}
    \end{pmatrix}
    \begin{pmatrix}
    \cos{\theta} & 0 & \sin{\theta} \\
    0 & 1 & 0\\
    -\sin{\theta} & 0 & \cos{\theta})
    \end{pmatrix}
    \begin{pmatrix}
    r\cos{\alpha}\\
    r\sin{\alpha}\\
    0
    \end{pmatrix}
    =
    \begin{pmatrix}
    r\cos{\alpha}\cos{\theta} \\
    r(\sin{\alpha}\cos{\phi} + \cos{\alpha}\sin{\theta}\sin{\phi}) \\
    r(\sin{\alpha}\sin{\phi} - \cos{\alpha}\sin{\theta}\cos{\phi})
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
    # Forces
    i_x = 5  # kgf
    i_y = 4  # kgf
    i_z = 0  # kgf
    m_ix = 0  # kgf-m
    m_iy = 0  # kgf-m
    m_iz = 0  # kgf-m

    # radius and angle (atm symmetrical)
    r = 0.8  # m
    angle_with_x = 45  # deg

    # coefficient for Thrust Momentum in motors
    c = 0.35  # "Adimensional"
    return angle_with_x, c, i_x, i_y, i_z, m_ix, m_iy, m_iz, r


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Process
    ### Transform SI units
    """
    )
    return


@app.cell
def _(G, angle_with_x, i_x, i_y, i_z, sp):
    # EF(x)
    F_x = i_x * 9.80665  # N
    # EF(y)
    F_y = i_y * 9.80665  # N
    # EF(z)
    F_z = G + i_z  # N

    # Radius angle with x
    alpha = dict()
    alpha[1] = sp.rad(angle_with_x)  # radians
    alpha[2] = sp.rad(180 + angle_with_x)  # radians
    alpha[3] = sp.rad(360 - angle_with_x)  # radians
    alpha[4] = sp.rad(90 + angle_with_x)  # radians
    return F_x, F_y, F_z, alpha


@app.cell
def _(mo):
    mo.md(r"""### Defining variables to solve""")
    return


@app.cell
def _(symbols):
    # Define symbolic variables
    T_1, T_2, T_3, T_4, phi, theta = symbols("T1 T2 T3 T4 phi theta")
    return T_1, T_2, T_3, T_4, phi, theta


@app.cell
def _(mo):
    mo.md(r"""### Force equations""")
    return


@app.cell
def _(Eq, F_x, F_y, F_z, T_1, T_2, T_3, T_4, cos, phi, sin, sqrt, theta):
    # Define equations
    equations = dict()

    equations[1] = Eq(
        T_1 * cos(theta) - T_2 * sin(theta) + T_3 * sin(theta) - T_4 * sin(theta),
        F_x,
    )
    equations[2] = Eq(
        -T_1 * sin(phi) + T_2 * sin(phi) + T_3 * sin(phi) - T_4 * sin(phi), F_y
    )
    equations[3] = Eq(
        (T_1 + T_2 + T_3 + T_4) * sqrt(1 - sin(theta) ** 2 - sin(phi) ** 2),
        F_z,
    )
    equations
    return (equations,)


@app.cell
def _(mo):
    mo.md(r"""### Arm rotated matrix""")
    return


@app.cell
def _(cos, np, sin, sp):
    def get_arm_rotated_matrix(r, alpha, theta, phi):
        r = sp.Matrix(np.array([[r * cos(alpha), r * sin(alpha), 0]]).T)
        return sp.rot_axis1(phi) @ sp.rot_ccw_axis2(theta) @ r
    return (get_arm_rotated_matrix,)


@app.cell
def _(alpha, get_arm_rotated_matrix, phi, r, theta):
    R = dict()
    R[1] = get_arm_rotated_matrix(r, alpha[1], theta, phi)
    R[2] = get_arm_rotated_matrix(r, alpha[2], theta, phi)
    R[3] = get_arm_rotated_matrix(r, alpha[3], theta, phi)
    R[4] = get_arm_rotated_matrix(r, alpha[4], theta, phi)
    R
    return (R,)


@app.cell
def _(mo):
    mo.md(r"""### Thrust vector for matrix equation""")
    return


@app.cell
def _(sin, sp, sqrt):
    def define_thrust_matrix(T, theta, phi):
        return sp.Matrix(
            [
                [T * sin(theta)],
                [T * sin(phi)],
                [T * sqrt(1 - sin(theta) ** 2 - sin(phi) ** 2)],
            ]
        )
    return (define_thrust_matrix,)


@app.cell
def _(T_1, T_2, T_3, T_4, define_thrust_matrix, phi, theta):
    T = dict()
    T[1] = define_thrust_matrix(T_1, theta, -phi)
    T[2] = define_thrust_matrix(T_2, -theta, phi)
    T[3] = define_thrust_matrix(T_3, theta, phi)
    T[4] = define_thrust_matrix(T_4, -theta, -phi)
    T
    return (T,)


@app.cell
def _(mo):
    mo.md(r"""### Calculate Momentum Matrix and add to equations""")
    return


@app.cell
def _(R, T, c, sp):
    m_matrix = sp.zeros(3, 1)

    for i in range(1, 4 + 1):
        temp_m = R[i].cross(T[i])
        print(temp_m)
        temp_m += c * T[i]
        m_matrix += temp_m
    m_matrix
    return (m_matrix,)


@app.cell
def _(Eq, equations, m_ix, m_iy, m_iz, m_matrix):
    equations[4] = Eq(m_matrix[0], m_ix)
    equations[5] = Eq(m_matrix[1], m_iy)
    equations[6] = Eq(m_matrix[2], m_iz)
    equations
    return


@app.cell
def _(mo):
    mo.md(r"""### Build system into A x = B format""")
    return


@app.cell
def _(sp):
    def build_symbolic_system(variables, equations):
        """
        Build A and B matrices symbolically

        variables: list of symbolic variables [x1, x2, ...]
        equations: list of Eq objects
        """
        n = len(variables)
        A = sp.zeros(n, n)  # n x n zero matrix
        B = sp.zeros(n, 1)  # n x 1 zero matrix

        for i, eq in enumerate(equations):
            # Extract coefficients for each variable
            for j, var in enumerate(variables):
                # Coefficient of variable in equation
                coeff = eq.lhs.coeff(var)
                A[i, j] = coeff

            # Constant term (move all variable terms to lhs, constant remains)
            B[i] = eq.rhs - (
                eq.lhs - sum(eq.lhs.coeff(var) * var for var in variables)
            )

        return A, B
    return (build_symbolic_system,)


@app.cell
def _(T_1, T_2, T_3, T_4, build_symbolic_system, equations, phi, theta):
    A, B = build_symbolic_system([T_1, T_2, T_3, T_4, phi, theta], list(equations.values()))
    A
    return A, B


@app.cell
def _(B):
    B
    return


@app.cell
def _(A, B, T_1, T_2, T_3, T_4, phi, solve, theta):
    solution = solve((A,B),T_1, T_2, T_3, T_4, phi, theta)
    solution
    return


if __name__ == "__main__":
    app.run()
