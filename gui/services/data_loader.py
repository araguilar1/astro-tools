"""
data_loader.py

Load and filter the periodic orbit CSV catalog.

Functions
---------
get_available_lagrange_points : List available Lagrange points
get_available_families        : List families for a given Lagrange point
get_available_modifiers       : List modifiers for a given LP + family
filter_orbits                 : Filter catalog rows by LP, family, modifier
orbit_row_to_ic               : Extract initial condition from a row
"""

import os
import csv

_CSV_PATH = os.path.join(
    os.path.dirname(__file__), "..", "..", "data", "Periodic Orbit Data", "periodicLagrangeOrbits.csv"
)

_catalog = None


def _load_catalog():
    global _catalog
    if _catalog is not None:
        return _catalog
    with open(_CSV_PATH, newline="") as f:
        reader = csv.DictReader(f)
        _catalog = list(reader)
    return _catalog


def get_available_lagrange_points():
    """Return sorted list of unique Lagrange point labels."""
    rows = _load_catalog()
    return sorted(set(r["LagrangePoint"] for r in rows))


def get_available_families(lagrange_point):
    """Return sorted list of unique family names for a Lagrange point."""
    rows = _load_catalog()
    return sorted(set(r["Family"] for r in rows if r["LagrangePoint"] == lagrange_point))


def get_available_modifiers(lagrange_point, family):
    """Return sorted list of unique modifiers for a LP + family combo."""
    rows = _load_catalog()
    mods = sorted(set(
        r["Modifier"] for r in rows
        if r["LagrangePoint"] == lagrange_point and r["Family"] == family
    ))
    return mods


def filter_orbits(lagrange_point, family, modifier=None):
    """
    Filter catalog by Lagrange point, family, and optionally modifier.

    Returns list of dicts with numeric fields converted to float.
    """
    rows = _load_catalog()
    filtered = []
    for r in rows:
        if r["LagrangePoint"] != lagrange_point:
            continue
        if r["Family"] != family:
            continue
        if modifier is not None and r["Modifier"] != modifier:
            continue
        filtered.append(r)
    return filtered


def orbit_row_to_ic(row):
    """
    Extract (state_6, period, jacobi, stability_index) from a catalog row.

    Parameters
    ----------
    row : dict
        A row from the catalog.

    Returns
    -------
    state_6 : list of float
        [x, y, z, vx, vy, vz]
    period : float
        Full orbital period T
    jacobi : float
        Jacobi constant value
    stability_index : float
        Stability index
    """
    state_6 = [
        float(row["x0"]),
        float(row["y0"]),
        float(row["z0"]),
        float(row["xDot0"]),
        float(row["yDot0"]),
        float(row["zDot0"]),
    ]
    period = float(row["T"])
    jacobi = float(row["JacobiValue"])
    stability_index = float(row["Stability_Index"])
    return state_6, period, jacobi, stability_index
