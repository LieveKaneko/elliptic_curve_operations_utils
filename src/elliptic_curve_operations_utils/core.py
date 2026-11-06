"""Core arithmetic for short Weierstrass elliptic curves over prime fields.

The curve is defined by the equation ``y^2 = x^3 + a*x + b (mod p)`` with
``p`` a prime.  Points are represented by two integers ``x`` and ``y`` in
``[0, p)``, or by the special value ``None`` (the point at infinity).

We deliberately keep the representation simple: a ``Point`` is a small
immutable value object, and a ``Curve`` holds the parameters and exposes
``add`` and ``mul``.  No attempt is made to use Jacobian coordinates or
constant-time arithmetic; those are optimisations for production use, not
for prototyping.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple


@dataclass(frozen=True)
class Point:
    """An affine point on a short Weierstrass curve.

    ``None`` is used to represent the point at infinity.  The class is
    frozen so that points can be used as dictionary keys and compared
    safely; the field values are stored as plain Python integers.
    """

    x: int
    y: int

    def __post_init__(self) -> None:
        if not isinstance(self.x, int) or not isinstance(self.y, int):
            raise TypeError("Point coordinates must be integers")

    def __repr__(self) -> str:  # pragma: no cover - trivial
        return f"Point({self.x}, {self.y})"


@dataclass(frozen=True)
class Curve:
    """A short Weierstrass curve ``y^2 = x^3 + a*x + b`` over ``F_p``.

    The discriminant ``4*a^3 + 27*b^2`` must be non-zero modulo ``p``;
    otherwise the curve is singular and the group law is not defined.
    """

    p: int
    a: int
    b: int

    def __post_init__(self) -> None:
        if self.p <= 1:
            raise ValueError("p must be a prime greater than 1")
        # Reduce a and b modulo p so callers can pass unreduced values.
        object.__setattr__(self, "a", self.a % self.p)
        object.__setattr__(self, "b", self.b % self.p)
        disc = (4 * pow(self.a, 3, self.p) + 27 * pow(self.b, 2, self.p)) % self.p
        if disc == 0:
            raise ValueError("singular curve: discriminant is zero modulo p")

    def is_on_curve(self, point: Optional[Point]) -> bool:
        """Return True if ``point`` satisfies the curve equation.

        The point at infinity (``None``) is considered on the curve by
        convention, since it is the identity element of the group.
        """
        if point is None:
            return True
        if not (0 <= point.x < self.p and 0 <= point.y < self.p):
            return False
        left = (point.y * point.y) % self.p
        right = (pow(point.x, 3, self.p) + self.a * point.x + self.b) % self.p
        return left == right

    def add(self, p1: Optional[Point], p2: Optional[Point]) -> Optional[Point]:
        """Add two points on the curve.

        Returns ``None`` for the point at infinity.  Raises ``ValueError``
        if either input is not on the curve; silently accepting invalid
        points would make debugging cryptographic prototypes harder.
        """
        if not self.is_on_curve(p1):
            raise ValueError(f"point {p1!r} is not on the curve")
        if not self.is_on_curve(p2):
            raise ValueError(f"point {p2!r} is not on the curve")

        if p1 is None:
            return p2
        if p2 is None:
            return p1

        if p1.x == p2.x and (p1.y + p2.y) % self.p == 0:
            # p1 == -p2, so the sum is the identity.
            return None

        if p1 == p2:
            # Point doubling: slope = (3*x^2 + a) / (2*y).
            # If y == 0 the point has order 2; the formula would divide by
            # zero, but in that case p1 == -p1 and the sum is the identity,
            # which was already handled above.
            num = (3 * p1.x * p1.x + self.a) % self.p
            den = (2 * p1.y) % self.p
        else:
            # Distinct points: slope = (y2 - y1) / (x2 - x1).
            num = (p2.y - p1.y) % self.p
            den = (p2.x - p1.x) % self.p

        inv = pow(den, -1, self.p)  # Python 3.8+ modular inverse.
        lam = (num * inv) % self.p
        x3 = (lam * lam - p1.x - p2.x) % self.p
        y3 = (lam * (p1.x - x3) - p1.y) % self.p
        return Point(x3, y3)

    def mul(self, k: int, point: Optional[Point]) -> Optional[Point]:
        """Return ``k * point`` using double-and-add.

        ``k`` may be zero or negative; negative scalars use the inverse of
        the point, which is well-defined because the group is abelian.
        """
        if not isinstance(k, int):
            raise TypeError("scalar must be an integer")
        if not self.is_on_curve(point):
            raise ValueError(f"point {point!r} is not on the curve")
        if k == 0 or point is None:
            return None
        if k < 0:
            neg = Point(point.x, (-point.y) % self.p)
            return self.mul(-k, neg)

        result: Optional[Point] = None
        addend: Optional[Point] = point
        while k:
            if k & 1:
                result = self.add(result, addend)
            addend = self.add(addend, addend)
            k >>= 1
        return result
