"""Elliptic curve point addition and scalar multiplication over prime fields.

This package provides a tiny, dependency-free implementation of the group law
for short Weierstrass curves ``y^2 = x^3 + a*x + b`` over a prime field.
It is intended for cryptographic prototyping, where you want to experiment
with curve parameters and point arithmetic without pulling in a large
elliptic-curve library.
"""

from .core import Curve, Point

__all__ = ["Curve", "Point"]
