import math

import pytest

from glyphipsi.angles import dihedral_degrees, normalize_degrees


def test_dihedral_degrees_known_right_angle():
    angle = dihedral_degrees((1, 0, 0), (0, 0, 0), (0, 1, 0), (0, 1, 1))

    assert angle == pytest.approx(-90.0)


@pytest.mark.parametrize(
    ("angle", "expected"),
    [
        (0.0, 0.0),
        (180.0, 180.0),
        (-180.0, -180.0),
        (181.0, -179.0),
        (-181.0, 179.0),
        (540.0, 180.0),
        (-540.0, -180.0),
    ],
)
def test_normalize_degrees(angle, expected):
    assert normalize_degrees(angle) == pytest.approx(expected)


def test_dihedral_degrees_rejects_degenerate_points():
    with pytest.raises(ValueError):
        dihedral_degrees((0, 0, 0), (0, 0, 0), (1, 0, 0), (2, 0, 0))


def test_dihedral_degrees_returns_finite_value():
    angle = dihedral_degrees((0, 0, 0), (1, 0, 0), (1, 1, 0), (2, 1, 1))

    assert math.isfinite(angle)
    assert -180.0 <= angle <= 180.0
