import numpy as np
import pytest

import Newton


def test_solve_scalar_quadratic():
    # F(x) = x^2 - 4, root at x = 2
    def eval_func(cdata, dx, X, f):
        f[0] = (X[0] + dx[0]) ** 2 - 4.0

    def tang_func(cdata, dx, X, J, counter):
        J[0, 0] = 2.0 * (X[0] + dx[0])

    root = Newton.Solve(1, np.array([3.0]), 1e-12, eval_func, tang_func, None)
    assert root[0] == pytest.approx(2.0, abs=1e-10)


def test_no_convergence_raises():
    def eval_func(cdata, dx, X, f):
        f[0] = (X[0] + dx[0]) ** 2 + 1.0  # no real root

    def tang_func(cdata, dx, X, J, counter):
        J[0, 0] = 2.0 * (X[0] + dx[0]) + 5.0  # keep J nonsingular

    with pytest.raises(Exception, match="no convergence"):
        Newton.Solve(1, np.array([1.0]), 1e-12, eval_func, tang_func, None,
                     maximum_iterations=10)
