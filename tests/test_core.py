import unittest

from elliptic_curve_operations_utils import Curve, Point


class TestCurveConstruction(unittest.TestCase):
    def test_rejects_singular_curve(self):
        # y^2 = x^3 + 3x + 2 over F_5 has discriminant 4*27 + 27*4 = 216 = 1 mod 5? No:
        # 4*3^3 + 27*2^2 = 108 + 108 = 216, 216 mod 5 = 1, so not singular.
        # Use a known singular curve: y^2 = x^3 + 1 over F_7 has 4*1 + 27*1 = 31 mod 7 = 3.
        # A singular example: a=0, b=0 gives discriminant 0.
        with self.assertRaises(ValueError):
            Curve(p=7, a=0, b=0)

    def test_rejects_non_prime_p(self):
        # We only check p > 1 here; primality is the caller's responsibility.
        with self.assertRaises(ValueError):
            Curve(p=1, a=1, b=1)

    def test_reduces_parameters_modulo_p(self):
        c = Curve(p=7, a=9, b=-1)
        self.assertEqual(c.a, 2)
        self.assertEqual(c.b, 6)


class TestPointValidation(unittest.TestCase):
    def setUp(self):
        # Curve y^2 = x^3 + 2x + 3 over F_97 (non-singular).
        self.curve = Curve(p=97, a=2, b=3)

    def test_identity_is_on_curve(self):
        self.assertTrue(self.curve.is_on_curve(None))

    def test_valid_point(self):
        # (3, 6) satisfies y^2 = 36, x^3+2x+3 = 27+6+3 = 36 mod 97.
        self.assertTrue(self.curve.is_on_curve(Point(3, 6)))

    def test_invalid_point(self):
        self.assertFalse(self.curve.is_on_curve(Point(3, 7)))

    def test_out_of_range_point(self):
        self.assertFalse(self.curve.is_on_curve(Point(97, 0)))


class TestAddition(unittest.TestCase):
    def setUp(self):
        self.curve = Curve(p=97, a=2, b=3)
        self.P = Point(3, 6)
        self.Q = Point(80, 10)
        # Verify Q is on the curve: 10^2 = 100 = 3 mod 97.
        # x^3+2x+3 = 512000 + 160 + 3 = 512163; 512163 mod 97 = 3? Let's check:
        # 97*5280 = 512160, remainder 3. Yes.

    def test_identity_left(self):
        self.assertEqual(self.curve.add(None, self.P), self.P)

    def test_identity_right(self):
        self.assertEqual(self.curve.add(self.P, None), self.P)

    def test_add_inverse(self):
        neg = Point(self.P.x, (-self.P.y) % self.curve.p)
        self.assertIsNone(self.curve.add(self.P, neg))

    def test_add_distinct(self):
        R = self.curve.add(self.P, self.Q)
        self.assertTrue(self.curve.is_on_curve(R))
        # Expected value computed by hand:
        # slope = (10-6)/(80-3) = 4/77 mod 97. 77^{-1} mod 97: 77*63=4851=0? 97*50=4850, so 77*63=1 mod 97.
        # slope = 4*63 = 252 = 58 mod 97.
        # x3 = 58^2 - 3 - 80 = 3364 - 83 = 3281 mod 97. 97*33=3201, remainder 80.
        # y3 = 58*(3 - 80) - 6 = 58*(-77) - 6 = -4466 - 6 = -4472 mod 97.
        # 97*46 = 4462, so -4472 = -10 = 87 mod 97.
        self.assertEqual(R, Point(80, 87))

    def test_double(self):
        R = self.curve.add(self.P, self.P)
        self.assertTrue(self.curve.is_on_curve(R))
        # slope = (3*9 + 2)/(12) = 29/12 mod 97. 12^{-1} = 89 (12*89=1068=1 mod 97).
        # slope = 29*89 = 2581 mod 97. 97*26=2522, remainder 59.
        # x3 = 59^2 - 6 = 3481 - 6 = 3475 mod 97. 97*35=3395, remainder 80.
        # y3 = 59*(3 - 80) - 6 = 59*(-77) - 6 = -4543 - 6 = -4549 mod 97.
        # 97*46 = 4462, -4549 + 4462 = -87 = 10 mod 97.
        self.assertEqual(R, Point(80, 10))

    def test_add_rejects_invalid_point(self):
        bad = Point(1, 1)
        with self.assertRaises(ValueError):
            self.curve.add(self.P, bad)


class TestScalarMultiplication(unittest.TestCase):
    def setUp(self):
        self.curve = Curve(p=97, a=2, b=3)
        self.P = Point(3, 6)

    def test_zero_scalar(self):
        self.assertIsNone(self.curve.mul(0, self.P))

    def test_one_scalar(self):
        self.assertEqual(self.curve.mul(1, self.P), self.P)

    def test_two_scalar(self):
        self.assertEqual(self.curve.mul(2, self.P), self.curve.add(self.P, self.P))

    def test_negative_scalar(self):
        neg_P = Point(self.P.x, (-self.P.y) % self.curve.p)
        self.assertEqual(self.curve.mul(-1, self.P), neg_P)
        self.assertEqual(self.curve.mul(-2, self.P), self.curve.mul(2, neg_P))

    def test_mul_rejects_invalid_point(self):
        with self.assertRaises(ValueError):
            self.curve.mul(2, Point(1, 1))

    def test_mul_rejects_non_integer(self):
        with self.assertRaises(TypeError):
            self.curve.mul(1.5, self.P)  # type: ignore[arg-type]


if __name__ == "__main__":
    unittest.main()
