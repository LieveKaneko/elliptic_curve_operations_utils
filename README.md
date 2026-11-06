# elliptic_curve_operations_utils

A dependency-free implementation of point addition and scalar multiplication on short Weierstrass curves over prime fields.

## Usage

```python
from elliptic_curve_operations_utils import Curve, Point

# y^2 = x^3 + 2x + 3 over F_97
curve = Curve(p=97, a=2, b=3)
P = Point(3, 6)
Q = Point(80, 10)

assert curve.is_on_curve(P)
R = curve.add(P, Q)
S = curve.mul(5, P)
```

## Why this exists

When prototyping cryptographic protocols, you often need point arithmetic on a small curve without pulling in a full library. This package gives you the group law directly, with explicit handling of the point at infinity and validation that points lie on the curve. The trade-off is speed and side-channel resistance: the implementation uses affine coordinates and Python integers, so it is not constant-time and not suitable for production secrets.

## Awkward edge

The point at infinity is represented by `None`, not by a special `Point` value. This keeps the type simple but means every call site must handle `None`; for example, `curve.mul(0, P)` returns `None`.

## API

- `Curve(p, a, b)` — constructs a curve. Raises `ValueError` if `p <= 1` or if the discriminant `4a^3 + 27b^2` is zero modulo `p`.
- `curve.is_on_curve(point)` — returns `True` if `point` is `None` or satisfies the curve equation.
- `curve.add(p1, p2)` — returns the sum, or `None` for the identity. Raises `ValueError` if either input is not on the curve.
- `curve.mul(k, point)` — returns `k * point` using double-and-add. `k` may be zero or negative. Raises `TypeError` if `k` is not an integer, and `ValueError` if `point` is not on the curve.
- `Point(x, y)` — an immutable affine point.
