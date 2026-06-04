"""Backbone dihedral angle math."""

from __future__ import annotations

import math
from collections.abc import Sequence

Point3D = Sequence[float]


def normalize_degrees(angle: float) -> float:
    """Normalize an angle to the documented -180 to +180 degree range."""
    normalized = ((angle + 180.0) % 360.0) - 180.0
    if math.isclose(normalized, -180.0, abs_tol=1e-12) and angle > 0:
        return 180.0
    if math.isclose(normalized, 0.0, abs_tol=1e-12):
        return 0.0
    return normalized


def dihedral_degrees(p0: Point3D, p1: Point3D, p2: Point3D, p3: Point3D) -> float:
    """Return the dihedral angle for four 3D points in degrees."""
    b0 = _subtract(p0, p1)
    b1 = _subtract(p2, p1)
    b2 = _subtract(p3, p2)

    b1_unit = _unit(b1)
    v = _subtract(b0, _scale(b1_unit, _dot(b0, b1_unit)))
    w = _subtract(b2, _scale(b1_unit, _dot(b2, b1_unit)))

    if _norm(v) == 0.0 or _norm(w) == 0.0:
        raise ValueError("Cannot calculate dihedral angle from degenerate points.")

    x = _dot(v, w)
    y = _dot(_cross(b1_unit, v), w)
    return normalize_degrees(math.degrees(math.atan2(y, x)))


def distance(p0: Point3D, p1: Point3D) -> float:
    """Return Euclidean distance between two 3D points."""
    return _norm(_subtract(p0, p1))


def _subtract(a: Point3D, b: Point3D) -> tuple[float, float, float]:
    return (float(a[0]) - float(b[0]), float(a[1]) - float(b[1]), float(a[2]) - float(b[2]))


def _scale(v: Point3D, factor: float) -> tuple[float, float, float]:
    return (float(v[0]) * factor, float(v[1]) * factor, float(v[2]) * factor)


def _dot(a: Point3D, b: Point3D) -> float:
    return float(a[0]) * float(b[0]) + float(a[1]) * float(b[1]) + float(a[2]) * float(b[2])


def _cross(a: Point3D, b: Point3D) -> tuple[float, float, float]:
    return (
        float(a[1]) * float(b[2]) - float(a[2]) * float(b[1]),
        float(a[2]) * float(b[0]) - float(a[0]) * float(b[2]),
        float(a[0]) * float(b[1]) - float(a[1]) * float(b[0]),
    )


def _norm(v: Point3D) -> float:
    return math.sqrt(_dot(v, v))


def _unit(v: Point3D) -> tuple[float, float, float]:
    norm = _norm(v)
    if norm == 0.0:
        raise ValueError("Cannot normalize a zero-length vector.")
    return (float(v[0]) / norm, float(v[1]) / norm, float(v[2]) / norm)
