"""Symmetric 3x3 tensor stored as a column (xx, yy, zz, xy, yz, zx).

Pure-Python replacement for the compiled ColTensor.pyd
(MeshTools/source/ColTensorModule.cpp).  The storage order matches the
``.sig`` stress files: sxx syy szz sxy syz sxz.
"""
import numpy as np

from . import Vec3D as _V


class ColTensor:
    __slots__ = ("ct",)
    __array_ufunc__ = None      # let numpy scalars defer to our operators

    def __init__(self, xx=0.0, yy=0.0, zz=0.0, xy=0.0, yz=0.0, zx=0.0):
        self.ct = [float(xx), float(yy), float(zz), float(xy), float(yz), float(zx)]

    # -- components ------------------------------------------------------
    def xx(self): return self.ct[0]
    def yy(self): return self.ct[1]
    def zz(self): return self.ct[2]
    def xy(self): return self.ct[3]
    def yz(self): return self.ct[4]
    def zx(self): return self.ct[5]

    def __repr__(self):
        return "%g %g %g %g %g %g" % tuple(self.ct)

    def matrix(self):
        xx, yy, zz, xy, yz, zx = self.ct
        return np.array([[xx, xy, zx], [xy, yy, yz], [zx, yz, zz]])

    # -- arithmetic ------------------------------------------------------
    def __add__(self, other):
        if not isinstance(other, ColTensor):
            raise TypeError("Both operands should be ColTensor instances")
        return ColTensor(*[a + b for a, b in zip(self.ct, other.ct)])

    def __sub__(self, other):
        if not isinstance(other, ColTensor):
            raise TypeError("Both operands should be ColTensor instances")
        return ColTensor(*[a - b for a, b in zip(self.ct, other.ct)])

    def __mul__(self, s):
        return ColTensor(*[a * s for a in self.ct])

    __rmul__ = __mul__

    def __truediv__(self, s):
        return ColTensor(*[a / s for a in self.ct])

    # -- scalar measures -------------------------------------------------
    def I1(self):
        return self.ct[0] + self.ct[1] + self.ct[2]

    def I2(self):
        xx, yy, zz, xy, yz, zx = self.ct
        return xx * yy + yy * zz + zz * xx - xy * xy - yz * yz - zx * zx

    def I3(self):
        return float(np.linalg.det(self.matrix()))

    def DeviatoricPart(self):
        m = self.I1() / 3.0
        xx, yy, zz, xy, yz, zx = self.ct
        return ColTensor(xx - m, yy - m, zz - m, xy, yz, zx)

    def J2(self):
        d = self.DeviatoricPart().ct
        return 0.5 * (d[0] ** 2 + d[1] ** 2 + d[2] ** 2) + d[3] ** 2 + d[4] ** 2 + d[5] ** 2

    def EffectiveStress(self):
        """von Mises equivalent stress."""
        return (3.0 * self.J2()) ** 0.5

    def Norm(self):
        xx, yy, zz, xy, yz, zx = self.ct
        return (xx * xx + yy * yy + zz * zz + 2.0 * (xy * xy + yz * yz + zx * zx)) ** 0.5

    def PrincipalValues(self):
        """Return ``(ColTensor(s1, s2, s3, 0, 0, 0), [e1, e2, e3])``.

        Principal values are sorted in DESCENDING order (s1 >= s2 >= s3) and
        ``e_i`` is the unit Vec3D eigenvector belonging to ``s_i``.  (Eigenvector
        signs are arbitrary, as with the original Jacobi routine.)
        """
        w, v = np.linalg.eigh(self.matrix())
        order = np.argsort(w)[::-1]
        w, v = w[order], v[:, order]
        evecs = [_V.Vec3D(v[0, i], v[1, i], v[2, i]) for i in range(3)]
        return ColTensor(w[0], w[1], w[2], 0.0, 0.0, 0.0), evecs
