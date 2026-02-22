"""
sidebar.py

Left panel layout: system selection, family browser, differential corrector inputs.
"""

import dash_bootstrap_components as dbc
from dash import html, dcc, dash_table


# All bodies from CelestialObject
BODY_OPTIONS = [
    "Sun", "Moon", "Mercury", "Venus", "Earth", "Mars",
    "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto",
    "Charon", "Nix", "Hydra", "Ganymede", "Titan",
    "Titania", "Ceres", "Enceladus", "Phobos", "Triton",
    "Callisto", "Europa",
]


def make_dropdown(id_, options=None, value=None, placeholder="Select...", clearable=False):
    return dcc.Dropdown(
        id=id_,
        options=[{"label": o, "value": o} for o in (options or [])],
        value=value,
        placeholder=placeholder,
        clearable=clearable,
        className="dash-dropdown",
    )


def create_sidebar():
    return dbc.Col(
        [
            # --- System Selection ---
            html.H5("System", className="sidebar-heading"),
            html.Label("Primary Body", className="input-label"),
            make_dropdown("primary-body", BODY_OPTIONS, "Earth"),
            html.Label("Secondary Body", className="input-label mt-2"),
            make_dropdown("secondary-body", BODY_OPTIONS, "Moon"),
            html.Div(id="system-info", className="system-info mt-2"),
            html.Hr(className="sidebar-hr"),

            # --- Family Selection ---
            html.H5("Orbit Family", className="sidebar-heading"),
            html.Label("Lagrange Point", className="input-label"),
            make_dropdown("lagrange-point-dropdown", ["L1", "L2", "L3", "L4", "L5"], "L1"),
            html.Label("Family", className="input-label mt-2"),
            make_dropdown("family-dropdown", placeholder="Select family..."),
            html.Label("Modifier", className="input-label mt-2"),
            make_dropdown("modifier-dropdown", placeholder="(none)", clearable=True),
            html.Div(id="orbit-count-badge", className="mt-2"),
            html.Label("Color By", className="input-label mt-2"),
            make_dropdown(
                "color-by-dropdown",
                ["Uniform", "Jacobi Constant", "Stability Index", "Period"],
                "Uniform",
            ),
            dbc.Button(
                "Render Family", id="render-family-btn", color="primary",
                className="mt-3 w-100",
            ),

            # --- Orbit Data Table ---
            html.Div(
                dash_table.DataTable(
                    id="orbit-data-table",
                    columns=[
                        {"name": "x0", "id": "x0", "type": "numeric"},
                        {"name": "z0", "id": "z0", "type": "numeric"},
                        {"name": "vy0", "id": "yDot0", "type": "numeric"},
                        {"name": "T", "id": "T", "type": "numeric"},
                        {"name": "Jacobi", "id": "JacobiValue", "type": "numeric"},
                        {"name": "SI", "id": "Stability_Index", "type": "numeric"},
                    ],
                    data=[],
                    row_selectable="single",
                    selected_rows=[],
                    style_table={"overflowX": "auto", "maxHeight": "200px", "overflowY": "auto"},
                    style_header={
                        "backgroundColor": "#0f0f2a",
                        "color": "#7090ff",
                        "fontWeight": "bold",
                        "border": "1px solid #1a1a3a",
                    },
                    style_cell={
                        "backgroundColor": "#0a0a1a",
                        "color": "#c0c0d0",
                        "border": "1px solid #1a1a3a",
                        "fontFamily": "JetBrains Mono, monospace",
                        "fontSize": "11px",
                        "padding": "4px 8px",
                    },
                    style_data_conditional=[
                        {"if": {"state": "selected"}, "backgroundColor": "#1a1a4a", "border": "1px solid #7090ff"},
                    ],
                ),
                className="mt-3",
            ),

            html.Hr(className="sidebar-hr"),

            # --- Differential Corrector ---
            html.H5("Differential Corrector", className="sidebar-heading"),
            *_corrector_inputs(),
            dbc.Button(
                "Correct Orbit", id="correct-orbit-btn", color="success",
                className="mt-3 w-100",
            ),
        ],
        width=3,
        className="sidebar",
    )


def _corrector_inputs():
    """Generate the 6 state inputs + period + fix-index inputs."""
    defaults = {
        "ic-x": 0.8234, "ic-y": 0.0, "ic-z": 0.0,
        "ic-vx": 0.0, "ic-vy": 0.126232, "ic-vz": 0.0,
    }
    state_labels = [
        ("x", "ic-x"), ("y", "ic-y"), ("z", "ic-z"),
        ("vx", "ic-vx"), ("vy", "ic-vy"), ("vz", "ic-vz"),
    ]
    elements = []
    row_items = []
    for label, id_ in state_labels:
        row_items.append(
            dbc.Col([
                html.Label(label, className="input-label-sm"),
                dbc.Input(
                    id=id_, type="number", value=defaults[id_],
                    className="ic-input", step="any",
                ),
            ], width=4)
        )
        if len(row_items) == 3:
            elements.append(dbc.Row(row_items, className="g-1 mb-1"))
            row_items = []

    elements.append(html.Label("Period (T)", className="input-label mt-2"))
    elements.append(dbc.Input(id="ic-period", type="number", value=2.743, className="ic-input", step="any"))

    elements.append(html.Label("Fix Front Indices", className="input-label mt-2"))
    elements.append(dbc.Input(id="fix-init-input", type="text", value="0,1,2,3,6", className="ic-input"))

    elements.append(html.Label("Fix Back Indices", className="input-label mt-2"))
    elements.append(dbc.Input(id="fix-end-input", type="text", value="1,3,5", className="ic-input"))

    return elements
