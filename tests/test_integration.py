import math

import pytest

from ddsim import Integration as I


def rhs(t, y, args):
    return y  # y' = y  ->  y(t) = exp(t)


def integrate(stepper, h, t_end=1.0):
    t, y = 0.0, 1.0
    while t < t_end - 1e-12:
        y = stepper(h, t, y, rhs, None)
        t += h
    return y


def test_euler_first_order():
    err = abs(integrate(I.Euler, 0.01) - math.e)
    assert 0.005 < err < 0.02


def test_rk4_fourth_order():
    e1 = abs(integrate(I.RK4, 0.1) - math.e)
    e2 = abs(integrate(I.RK4, 0.05) - math.e)
    assert e1 < 1e-5
    assert 12 < e1 / e2 < 20  # halving h cuts the error by ~2**4


def test_midpoint_second_order():
    e1 = abs(integrate(I.MidPoint, 0.1) - math.e)
    e2 = abs(integrate(I.MidPoint, 0.05) - math.e)
    assert 3 < e1 / e2 < 5


def test_rk5_ck_accurate_and_returns_next_step():
    y, step_taken, next_step = I.RK5_CK(0.1, 0.0, 1.0, rhs, None, 1e-8)
    assert y == pytest.approx(math.exp(step_taken), rel=1e-8)
    assert next_step > 0


def test_vector_helpers():
    assert I.plus([1, 2], [3, 4]) == [4, 6]
    assert I.minus([3, 4], [1, 2]) == [2, 2]
    assert I.star([1, 2], 3) == [3, 6]
    assert I.norm([3, 4]) == pytest.approx(5.0)
    assert I.maxnorm([3, -4]) == pytest.approx(4.0)
