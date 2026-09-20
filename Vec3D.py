"""Three-element vector.

Pure-Python replacement for the compiled Vec3D.pyd (MeshTools/source/Vec3DModule.cpp).
Semantics follow the C++ module:

* ``v + v``, ``v - v``; ``v + s``, ``v - s``, ``s - v`` apply the scalar to every component
* ``v * s`` / ``s * v`` scale; ``v * v`` is the DOT product (a float)
* ``v / s`` divides by a number
* ``Normalize()`` returns a NEW unit vector; ``Max`` / ``Min`` mutate in place
* ``repr`` is ``"%g %g %g"``
"""
import math
from numbers import Real


def _is_num(o):
    return isinstance(o, Real)


class Vec3D:
    __slots__ = ("_v",)
    # Make numpy scalars/arrays defer to our operators: without this,
    # ``np.float64(2) * v`` sees a sequence and returns an ndarray.
    __array_ufunc__ = None

    def __init__(self, x, y, z):
        self._v = [float(x), float(y), float(z)]

    # -- accessors -------------------------------------------------------
    def x(self):
        return self._v[0]

    def y(self):
        return self._v[1]

    def z(self):
        return self._v[2]

    def __getitem__(self, i):
        return self._v[i]

    def __setitem__(self, i, value):
        self._v[i] = float(value)

    def __len__(self):
        return 3

    def __iter__(self):
        return iter(self._v)

    def __repr__(self):
        return "%g %g %g" % tuple(self._v)

    __str__ = __repr__

    # -- arithmetic ------------------------------------------------------
    def __add__(self, other):
        a = self._v
        if isinstance(other, Vec3D):
            b = other._v
            return Vec3D(a[0] + b[0], a[1] + b[1], a[2] + b[2])
        if _is_num(other):
            return Vec3D(a[0] + other, a[1] + other, a[2] + other)
        return NotImplemented

    __radd__ = __add__

    def __sub__(self, other):
        a = self._v
        if isinstance(other, Vec3D):
            b = other._v
            return Vec3D(a[0] - b[0], a[1] - b[1], a[2] - b[2])
        if _is_num(other):
            return Vec3D(a[0] - other, a[1] - other, a[2] - other)
        return NotImplemented

    def __rsub__(self, other):
        if _is_num(other):
            a = self._v
            return Vec3D(other - a[0], other - a[1], other - a[2])
        return NotImplemented

    def __mul__(self, other):
        a = self._v
        if isinstance(other, Vec3D):
            b = other._v
            return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
        if _is_num(other):
            return Vec3D(a[0] * other, a[1] * other, a[2] * other)
        return NotImplemented

    def __rmul__(self, other):
        if _is_num(other):
            a = self._v
            return Vec3D(other * a[0], other * a[1], other * a[2])
        return NotImplemented

    def __truediv__(self, other):
        if _is_num(other):
            a = self._v
            return Vec3D(a[0] / other, a[1] / other, a[2] / other)
        return NotImplemented

    def __neg__(self):
        a = self._v
        return Vec3D(-a[0], -a[1], -a[2])

    def __pos__(self):
        return Vec3D(*self._v)

    # in-place operators mutate, as in the C++ module
    def __iadd__(self, other):
        r = self + other
        if r is NotImplemented:
            return r
        self._v[:] = r._v
        return self

    def __isub__(self, other):
        r = self - other
        if r is NotImplemented:
            return r
        self._v[:] = r._v
        return self

    def __imul__(self, other):
        if not _is_num(other):
            raise TypeError("Can only multiply a Vec3D by a number")
        self._v[:] = [c * other for c in self._v]
        return self

    def __itruediv__(self, other):
        if not _is_num(other):
            raise TypeError("Can only divide a Vec3D by a number")
        self._v[:] = [c / other for c in self._v]
        return self

    # -- methods ---------------------------------------------------------
    def Magnitude(self):
        a = self._v
        return math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])

    def Normalize(self):
        """Return a new unit vector (does not modify self)."""
        length = self.Magnitude()
        if length == 0.0:
            raise TypeError("Cannot normalize a zero length Vec3D")
        a = self._v
        return Vec3D(a[0] / length, a[1] / length, a[2] / length)

    def Max(self, other):
        """Component-wise max, stored in self."""
        for i in range(3):
            if other._v[i] > self._v[i]:
                self._v[i] = other._v[i]
        return self

    def Min(self, other):
        """Component-wise min, stored in self."""
        for i in range(3):
            if other._v[i] < self._v[i]:
                self._v[i] = other._v[i]
        return self


def CrossProd(a, b):
    a, b = a._v, b._v
    return Vec3D(a[1] * b[2] - a[2] * b[1],
                 a[2] * b[0] - a[0] * b[2],
                 a[0] * b[1] - a[1] * b[0])


def TripleProd(a, b, c):
    """(a x b) . c"""
    return CrossProd(a, b) * c
