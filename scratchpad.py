import marimo

__generated_with = "0.17.0"
app = marimo.App(
    width="medium",
    auto_download=["ipynb"],
    sql_output="lazy-polars",
)


@app.cell
def _():
    return


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
    return


@app.cell
def _(A, Eq, linsolve, sp, symbols):
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
    solution = linsolve((A_1, B_1), x1, x2, x3)
    print(f"\nSolution: {solution}")
    return (B_1,)


if __name__ == "__main__":
    app.run()
