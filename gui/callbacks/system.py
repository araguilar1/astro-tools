"""
system.py

Callbacks for body pair selection and system parameter computation.
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from dash import Input, Output, callback, html
from astro import CelestialObject, characteristicQuantities, lagrangePoints


@callback(
    Output("system-params", "data"),
    Output("system-info", "children"),
    Input("primary-body", "value"),
    Input("secondary-body", "value"),
)
def update_system(primary_name, secondary_name):
    if not primary_name or not secondary_name or primary_name == secondary_name:
        return {}, html.Div("Select two different bodies.", className="text-warning")

    try:
        p1 = CelestialObject(primary_name)
        p2 = CelestialObject(secondary_name)
        cq = characteristicQuantities(p1, p2)
    except Exception as e:
        return {}, html.Div(f"Error: {e}", className="text-danger")

    mu = cq["mu"]
    lstar = cq["lstar"]
    tstar = cq["tstar"]

    lps = lagrangePoints(mu)

    params = {
        "mu": mu,
        "lstar": lstar,
        "tstar": tstar,
        "lagrange_points": lps.tolist(),
    }

    info = html.Div([
        html.Div(f"\u03bc = {mu:.10f}", className="sys-param"),
        html.Div(f"l* = {lstar:.2f} km", className="sys-param"),
        html.Div(f"t* = {tstar:.2f} s", className="sys-param"),
    ])

    return params, info
