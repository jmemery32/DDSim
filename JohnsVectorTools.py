"""Element-wise arithmetic on equal-length sequences; results are lists.

Pure-Python replacement for the compiled JohnsVectorTools.pyd (source was not
recovered; behavior inferred from its call sites in DamClass, DamMo and PDF).
"""


def Plus(a, b):
    return [x + y for x, y in zip(a, b)]


def Minus(a, b):
    return [x - y for x, y in zip(a, b)]


def Star(a, b):
    """Element-wise product."""
    return [x * y for x, y in zip(a, b)]


def Divide(a, b):
    """Element-wise quotient a/b."""
    return [x / y for x, y in zip(a, b)]


def ScalarMult(a, s):
    return [x * s for x in a]
