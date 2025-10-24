import marimo

__generated_with = "0.17.0"
app = marimo.App(
    width="medium",
    auto_download=["ipynb"],
    sql_output="lazy-polars",
)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# Mathematical force analysis of VTOL""")
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _():
    from math import radians
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
    return Eq, linsolve, sin, sp, sqrt, symbols


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Starting point""")
    return


@app.cell(hide_code=True)
def _(mo):
    # Gravitational force
    m = mo.ui.slider(
        10,
        40,
        0.1,
        25.5,
        label="UAV Mass (kg):",
        debounce=True,
        show_value=True,
        full_width=True,
    )  # kg o kgf
    g = 9.81  # m/s^2

    mo.hstack([m])
    return g, m


@app.cell
def _(g, m):
    G = m.value * g
    return (G,)


@app.cell
def _(G, g, np, radians):
    # Pitch
    th = 10
    # Roll
    ph = 0

    # First iteration
    if abs(th) > 90 or abs(ph) > 90:
        raise AttributeError("Any angle can be over 90 degrees")

    firstT = G / np.sqrt(1 - np.sin(radians(th)) ** 2 - np.sin(radians(ph)) ** 2)
    print(f"The thrust per motor should be {firstT / (g * 4):.3f} kgf")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Base concepts


    ### Momentum

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
    ### Rotation Matrix
    With clockwise rotation around X-axis (roll $\phi$) and clockwise rotation around Y-axis (pitch $\theta$) with a position vector $\vec{r}$ with angle  $\alpha$ : $\vec{R} = \Phi\times\Theta\times\vec{r}$. It differs with us for one direction.

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
    mo.md(
        r"""
    ### Thrust vector 

    Usually orthogonal to the XY plane, upwards, with rotation given two euler angle, Roll = $\phi$ and Pitch = $\theta$.

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


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Main  Equations

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
    \vec{M_i} + \vec{r_i}\times\vec{i} &= c_{ccw}\vec{T_1} + \vec{r_1}\times\vec{T_1} \\
              &+ c_{ccw}\vec{T_2} + \vec{r_2}\times\vec{T_2} \\
              &+ c_{cw}\vec{T_3} + \vec{r_3}\times\vec{T_3} \\
              &+ c_{cw}\vec{T_4} + \vec{r_4}\times\vec{T_4}
    \end{align}
    $$
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Inputs""")
    return


@app.cell(hide_code=True)
def _(mo):
    # Angles
    roll = mo.ui.slider(
        start=-90,
        stop=90,
        step=0.05,
        value=0,
        label="Roll:",
        debounce=True,
        show_value=True,
        full_width=True,
    )  # degrees
    pitch = mo.ui.slider(
        start=-90,
        stop=90,
        step=0.1,
        value=0,
        label="Pitch:",
        debounce=True,
        show_value=True,
        full_width=True,
    )  # degrees

    return pitch, roll


@app.cell(hide_code=True)
def _(mo, pitch, roll):
    mo.vstack([roll, pitch])
    return


@app.cell(hide_code=True)
def _(mo):
    # Forces
    switch = mo.ui.switch(value=True, label="kg")
    i_x = mo.ui.number(value=0, step=0.01, debounce=True)  # kgf
    i_y = mo.ui.number(value=0, step=0.01, debounce=True)  # kgf
    i_z = mo.ui.number(value=0, step=0.01, debounce=True)  # kgf
    # Momentum
    m_switch = mo.ui.switch(value=True, label="kg-m")
    i_m_x = mo.ui.number(value=0, step=0.01, debounce=True)  # kgf-m
    i_m_y = mo.ui.number(value=0, step=0.01, debounce=True)  # kgf-m
    i_m_z = mo.ui.number(value=0, step=0.01, debounce=True)  # kgf-m

    # motors coefficient for Thrust Momentum
    m_coeff = mo.ui.slider(
        start=0.05,
        stop=0.85,
        step=0.05,
        value=0.05,
        debounce=True,
        show_value=True,
        full_width=True,
    )  # "Adimensional"
    return i_m_x, i_m_y, i_m_z, i_x, i_y, i_z, m_coeff, m_switch, switch


@app.cell(hide_code=True)
def _(i_m_x, i_m_y, i_m_z, i_x, i_y, i_z, m_coeff, m_switch, mo, switch):
    mo.accordion(
        {
            "X force input": mo.hstack([i_x, switch]),
            "Y force input": mo.hstack([i_y, switch]),
            "Z force input": mo.hstack([i_z, switch]),
            "X momentum input": mo.hstack([i_m_x, m_switch]),
            "Y momentum input": mo.hstack([i_m_y, m_switch]),
            "Z momentum input": mo.hstack([i_m_z, m_switch]),
            "Motors Momentum coefficient": m_coeff,
        },
        multiple=True,
    )
    return


@app.cell(hide_code=True)
def _(mo):
    # Forces applied location
    r_x_i = mo.ui.slider(
        -2, 2, 0.1, 0, label="F(x)", show_value=True, full_width=True
    )
    r_y_i = mo.ui.slider(
        -3.6, 3.6, 0.1, 0, label="F(y)", show_value=True, full_width=True
    )
    r_z_i = mo.ui.slider(
        -1, 1, 0.1, 0, label="F(z)", show_value=True, full_width=True
    )
    # Gravity applied location
    r_x_g = mo.ui.slider(
        -2, 2, 0.1, 0, label="G(x)", show_value=True, full_width=True
    )
    r_y_g = mo.ui.slider(
        -3.6, 3.6, 0.1, 0, label="G(y)", show_value=True, full_width=True
    )
    r_z_g = mo.ui.slider(
        -1, 1, 0.1, 0, label="G(z)", show_value=True, full_width=True
    )
    return r_x_g, r_x_i, r_y_g, r_y_i, r_z_g, r_z_i


@app.cell(hide_code=True)
def _(mo, r_x_g, r_x_i, r_y_g, r_y_i, r_z_g, r_z_i):
    mo.vstack(
        [
            mo.md("Force location"),
            r_x_i,
            r_y_i,
            r_z_i,
            mo.md("Gravity location"),
            r_x_g,
            r_y_g,
            r_z_g,
        ],
        align="center",
    )
    return


@app.cell(hide_code=True)
def _(mo):
    # VTOL motors locations
    # r = 0.8  # m
    # angle_with_x = 45  # deg
    r_x_1 = mo.ui.slider(
        0.1, 1, 0.05, 0.4, orientation="vertical", show_value=True
    )
    r_y_1 = mo.ui.slider(0.1, 1, 0.05, 0.4, debounce=True, show_value=True)
    r_z_1 = mo.ui.slider(
        -0.5, 0.5, 0.05, 0, label="z1", orientation="vertical", show_value=True
    )

    r_x_2 = mo.ui.slider(
        -1, -0.1, 0.05, -0.4, orientation="vertical", show_value=True
    )
    r_y_2 = mo.ui.slider(-1, -0.1, 0.05, -0.4, debounce=True, show_value=True)
    r_z_2 = mo.ui.slider(
        -0.5, 0.5, 0.05, 0, label="z2", orientation="vertical", show_value=True
    )

    r_x_3 = mo.ui.slider(
        0.1, 1, 0.05, 0.4, orientation="vertical", show_value=True
    )
    r_y_3 = mo.ui.slider(-1, -0.1, 0.05, -0.4, debounce=True, show_value=True)
    r_z_3 = mo.ui.slider(
        -0.5, 0.5, 0.05, 0, label="z3", orientation="vertical", show_value=True
    )

    r_x_4 = mo.ui.slider(
        -1, -0.1, 0.05, -0.4, orientation="vertical", show_value=True
    )
    r_y_4 = mo.ui.slider(0.1, 1, 0.05, 0.4, debounce=True, show_value=True)
    r_z_4 = mo.ui.slider(
        -0.5, 0.5, 0.05, 0, label="z4", orientation="vertical", show_value=True
    )
    return (
        r_x_1,
        r_x_2,
        r_x_3,
        r_x_4,
        r_y_1,
        r_y_2,
        r_y_3,
        r_y_4,
        r_z_1,
        r_z_2,
        r_z_3,
        r_z_4,
    )


@app.cell(hide_code=True)
def _(
    mo,
    r_x_1,
    r_x_2,
    r_x_3,
    r_x_4,
    r_y_1,
    r_y_2,
    r_y_3,
    r_y_4,
    r_z_1,
    r_z_2,
    r_z_3,
    r_z_4,
):
    mo.hstack(
        [
            mo.vstack(
                [mo.md("^ X T3"), r_x_3, r_y_3, mo.md("T2"), r_y_2, r_x_2],
                align="end",
            ),
            mo.vstack(
                [mo.md("T1"), r_x_1, r_y_1, mo.md("T4 Y >"), r_y_4, r_x_4],
                align="start",
            ),
            r_z_1,
            r_z_2,
            r_z_3,
            r_z_4,
        ]
    )
    return


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
def _():
    N_kg = 9.80665
    return (N_kg,)


@app.cell
def _(N_kg, i_x, i_y, i_z, switch):
    # F(x)
    F_x = i_x.value * N_kg if switch.value else i_x.value  # N
    # F(y)
    F_y = i_y.value * N_kg if switch.value else i_y.value  # N
    # F(z)
    F_z = i_z.value * N_kg if switch.value else i_z.value  # N
    return F_x, F_y, F_z


@app.cell
def _(N_kg, i_m_x, i_m_y, i_m_z, m_switch):
    # M(x)
    m_x = i_m_x.value * N_kg if m_switch.value else i_m_x.value  # Nm
    # M(y)
    m_y = i_m_y.value * N_kg if m_switch.value else i_m_y.value  # Nm
    # M(z)
    m_z = i_m_z.value * N_kg if m_switch.value else i_m_z.value  # Nm
    # m_z = i_m_x.value * m_switch.value # (1 for Nm, 9.80665 for kgf, etc)
    return m_x, m_y, m_z


@app.cell
def _(m_coeff):
    m_c_cw = -m_coeff.value
    m_c_ccw = m_coeff.value
    return m_c_ccw, m_c_cw


@app.cell
def _(pitch, roll, sp):
    # Radius angle with x
    # alpha = dict()
    # alpha[1] = sp.rad(angle_with_x)  # radians
    # alpha[2] = sp.rad(180 + angle_with_x)  # radians
    # alpha[3] = sp.rad(360 - angle_with_x)  # radians
    # alpha[4] = sp.rad(90 + angle_with_x)  # radians

    # Euler Angles
    phi = sp.rad(roll.value)  # Roll in radians
    theta = sp.rad(pitch.value)  # Pitch in radians
    return phi, theta


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Calculate inputs""")
    return


@app.cell
def _(F_z, G):
    # EF(z)
    EF_z = G + F_z
    return (EF_z,)


@app.cell
def _(
    F_x,
    F_y,
    F_z,
    G,
    m_x,
    m_y,
    m_z,
    np,
    r_x_g,
    r_x_i,
    r_y_g,
    r_y_i,
    r_z_g,
    r_z_i,
    sp,
):
    i_r = sp.Matrix(np.array([[r_x_i.value, r_y_i.value, r_z_i.value]]).T)
    m_i = sp.Matrix(np.array([[F_x, F_y, F_z]]).T)
    M_i = i_r.cross(m_i)
    # M_g
    i_g = sp.Matrix(np.array([[r_x_g.value, r_y_g.value, r_z_g.value]]).T)
    m_g = sp.Matrix(np.array([[0, 0, G]]).T)
    M_g = i_g.cross(m_g)
    # EM
    EM = sp.Matrix(np.array([[m_x, m_y, m_z]]).T) + M_g + M_i
    M_x = EM[0]
    M_y = EM[1]
    M_z = EM[2]
    EM
    return M_x, M_y, M_z


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Defining variables to solve""")
    return


@app.cell
def _(symbols):
    # Define symbolic variables
    T_1, T_2, T_3, T_4 = symbols("T1 T2 T3 T4")
    return T_1, T_2, T_3, T_4


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Force equations""")
    return


@app.cell(hide_code=True)
def _(EF_z, Eq, F_x, F_y, T_1, T_2, T_3, T_4, phi, sin, sqrt, theta):
    # Define equations
    equations = dict()

    equations[1] = Eq(
        T_1 * sin(theta) - T_2 * sin(theta) + T_3 * sin(theta) - T_4 * sin(theta),
        F_x,
    )
    equations[2] = Eq(
        -T_1 * sin(phi) + T_2 * sin(phi) + T_3 * sin(phi) - T_4 * sin(phi), F_y
    )
    equations[3] = Eq(
        (T_1 + T_2 + T_3 + T_4) * sqrt(1 - sin(theta) ** 2 - sin(phi) ** 2),
        EF_z,
    )
    equations
    return (equations,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### VTOL location matrices with rotations""")
    return


@app.cell(hide_code=True)
def _(np, sp):
    def get_arm_rotated_matrix(r_x, r_y, r_z, theta, phi):
        r = sp.Matrix(np.array([[r_x, r_y, r_z]]).T)
        return sp.rot_axis1(phi) @ sp.rot_ccw_axis2(theta) @ r
    return (get_arm_rotated_matrix,)


@app.cell(hide_code=True)
def _(
    get_arm_rotated_matrix,
    phi,
    r_x_1,
    r_x_2,
    r_x_3,
    r_x_4,
    r_y_1,
    r_y_2,
    r_y_3,
    r_y_4,
    r_z_1,
    r_z_2,
    r_z_3,
    r_z_4,
    theta,
):
    R = dict()
    R[1] = get_arm_rotated_matrix(
        r_x_1.value, r_y_1.value, r_z_1.value, theta, phi
    )
    R[2] = get_arm_rotated_matrix(
        r_x_2.value, r_y_2.value, r_z_2.value, theta, phi
    )
    R[3] = get_arm_rotated_matrix(
        r_x_3.value, r_y_3.value, r_z_3.value, theta, phi
    )
    R[4] = get_arm_rotated_matrix(
        r_x_4.value, r_y_4.value, r_z_4.value, theta, phi
    )
    R
    return (R,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Thrust vector for matrix equation""")
    return


@app.cell(hide_code=True)
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


@app.cell(hide_code=True)
def _(T_1, T_2, T_3, T_4, define_thrust_matrix, phi, theta):
    T = dict()
    T[1] = define_thrust_matrix(T_1, theta, -phi)
    T[2] = define_thrust_matrix(T_2, -theta, phi)
    T[3] = define_thrust_matrix(T_3, theta, phi)
    T[4] = define_thrust_matrix(T_4, -theta, -phi)
    T
    return (T,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Calculate Momentum Matrix and add to equations""")
    return


@app.cell(hide_code=True)
def _(R, T, m_c_ccw, m_c_cw, sp):
    m_matrix = sp.zeros(3, 1)

    for i in range(1, 4 + 1):
        temp_m = R[i].cross(T[i])
        print(temp_m)
        if i <= 2:
            temp_m += m_c_ccw * T[i]
        else:
            temp_m += m_c_cw * T[i]
        m_matrix += temp_m
    m_matrix
    return (m_matrix,)


@app.cell(hide_code=True)
def _(Eq, M_x, M_y, M_z, equations, m_matrix, sp):
    equations[4] = Eq(m_matrix[0], M_x)
    equations[5] = Eq(m_matrix[1], M_y)
    equations[6] = Eq(m_matrix[2], M_z)

    parsed_equations = list()

    for k in equations:
        if equations.get(k) not in (True, sp.true, False, sp.false):
            parsed_equations.append(equations[k])
    parsed_equations
    return (parsed_equations,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Outputs
    ### Build system into $A\vec{x} = B$ format
    """
    )
    return


@app.cell(hide_code=True)
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
def _(T_1, T_2, T_3, T_4, build_symbolic_system, parsed_equations):
    if len(parsed_equations) > 4:
        A, B = build_symbolic_system([T_1, T_2, T_3, T_4], parsed_equations[1:5])
    else:
        A, B = build_symbolic_system([T_1, T_2, T_3, T_4], parsed_equations)
    A
    return A, B


@app.cell
def _(B):
    B
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Symbolic resolution""")
    return


@app.cell
def _(A, B, T_1, T_2, T_3, T_4, linsolve):
    solution = linsolve((A, B), T_1, T_2, T_3, T_4)
    solution
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Numeric resolution""")
    return


@app.cell
def _(A, B, np):
    A_numeric = np.array(A, dtype=np.float64)
    B_numeric = np.array(B, dtype=np.float64)
    print(B)
    A
    return A_numeric, B_numeric


@app.cell
def _(A_numeric, B_numeric, N_kg, np):
    numeric_solution = np.linalg.solve(A_numeric, B_numeric).squeeze()
    for n in range(len(numeric_solution)):
        print(
            f"T{n + 1} has a magnitude of {numeric_solution[n] / N_kg:.3f}kg."
        )
    return


if __name__ == "__main__":
    app.run()
