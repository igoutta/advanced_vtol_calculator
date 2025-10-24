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
    from math import sin, cos, radians, degrees, pi
    return degrees, pi, radians, sin


@app.cell
def _(mo):
    mo.md(r"""## Degrees Sine Test""")
    return


@app.cell
def _(radians, sin):
    for iota in range(30):
        print(
            f"Comparison of {iota}° in radians: {radians(iota):.3f} and {sin(radians(iota)):.3f}"
        )
    return


@app.cell
def _(mo):
    mo.md(r"""## Numpy matrices math test""")
    return


@app.cell
def _():
    import numpy as np
    return (np,)


@app.cell
def _(np):
    A = np.array([[9, 8, 7]])
    A.shape
    return (A,)


@app.cell
def _(A):
    A.size
    return


@app.cell
def _(np):
    B = np.array([[3, 4, 5]]).T
    B
    return (B,)


@app.cell
def _(A, B):
    A @ B
    return


@app.cell
def _(A, B, np):
    Momento = np.array([np.ones(A.size), A.squeeze(), B.squeeze()])
    Momento
    return


@app.cell
def _(A, B, np):
    np.linalg.cross(A.squeeze(), B.squeeze())
    return


@app.cell
def _(B, np):
    2 * np.eye(3) @ B
    return


@app.cell
def _(mo):
    mo.md(r"""## Sympy symbolic equations test""")
    return


@app.cell
def _():
    import sympy as sp
    from sympy import Eq, linsolve, solve, symbols
    return Eq, linsolve, solve, sp, symbols


@app.cell
def _(Eq, solve, symbols):
    # Define symbolic variables
    x, y, z = symbols("x y z")

    # Create equations
    eq1 = Eq(2 * x + 3 * y - z, 5)
    eq2 = Eq(x - y + 2 * z, 3)
    eq3 = Eq(3 * x + y - z, 7)

    print("Equations:")
    print(f"1: {eq1}")
    print(f"2: {eq2}")
    print(f"3: {eq3}")

    # Solve symbolically
    solution = solve([eq1, eq2, eq3], [x, y, z])
    print(f"\nSymbolic solution: {solution}")
    return x, y, z


@app.cell
def _(Eq, linsolve, sp, symbols):
    # Define variables
    x1, x2, x3 = symbols("x1 x2 x3")

    # Coefficient matrix with symbols
    A_1 = sp.Matrix([[2, 1, -1], [-3, -1, 2], [-2, 1, 2]])

    # Variable vector
    X = sp.Matrix([x1, x2, x3])

    # Constant vector
    B_1 = sp.Matrix([8, -11, -3])

    # Matrix equation: A * X = B
    matrix_eq = Eq(A_1 * X, B_1)
    print("Matrix equation:")
    print(matrix_eq)

    # Solve
    solution_m = linsolve((A_1, B_1), x1, x2, x3)
    print(f"\nSolution: {solution_m}")
    return


@app.cell
def _(Eq, solve, symbols, x, y, z):
    # Define variables and parameters
    # x, y, z = symbols('x y z')
    m, n, p = symbols('m n p')  # Parameters

    # Create equations with parameters
    eq_1 = Eq(m*x + n*y - z, 5)
    eq_2 = Eq(2*x - y + p*z, 3)
    eq_3 = Eq(x + 2*y - 3*z, 7)

    print("Parameterized equations:")
    print(f"1: {eq_1}")
    print(f"2: {eq_2}")
    print(f"3: {eq_3}")

    # Solve in terms of parameters
    solution_p = solve([eq_1, eq_2, eq_3], [x, y, z])
    print("\nSolution in terms of parameters:")
    for var, expr in solution_p.items():
        print(f"{var} = {expr}")
    return


@app.cell
def _(Eq, linsolve, sp, x, y, z):
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
            B[i] = eq.rhs - (eq.lhs - sum(eq.lhs.coeff(var)*var for var in variables))

        return A, B

    # Example usage
    equations = [
        Eq(2*x + 3*y - z, 5),
        Eq(x - y + 2*z, 3),
        Eq(3*x + y - z, 7)
    ]

    A_2, B_2 = build_symbolic_system([x, y, z], equations)

    print("Coefficient matrix A:")
    print(A_2)
    print("\nConstant vector B:")
    print(B_2)

    # Solve using matrix form
    solution_m2 = linsolve((A_2, B_2), x, y, z)
    print(f"\nSolution: {solution_m2}")
    return (A_2,)


@app.cell
def _(pi):
    pi
    return


@app.cell
def _(degrees, pi, sp):
    degrees(sp.pi) == degrees(pi)
    return


@app.cell
def _(A_2):
    len(A_2)
    return


@app.cell
def _(mo):
    mo.md(r"""Rotations""")
    return


@app.cell
def _(np, sp, symbols):
    phi, theta = symbols('phi theta')
    r = sp.Matrix(np.array([[0.8*sp.cos(sp.pi/4), 0.8*sp.sin(sp.pi/4),0]]).T)
    r_1 = sp.rot_axis1(phi)@sp.rot_ccw_axis2(theta)@r
    r_1
    return


@app.cell
def _():
    from sympy import false

    print(f"This is: {false == False}")

    alpha = {
        1: True,
        2: false,
        3: False,
    }

    alist = list()
    for k in alpha:
        if alpha.get(k) is not False:
            alist.append(k)
    alist
    return


if __name__ == "__main__":
    app.run()
