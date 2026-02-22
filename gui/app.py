"""
app.py

CR3BP Explorer — Plotly Dash web application entry point.

Usage
-----
    python gui/app.py
    # Open http://127.0.0.1:8050 in browser
"""

import sys
import os

# Ensure project root is on the path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))

import dash
import dash_bootstrap_components as dbc
from dash import html, dcc

from components.sidebar import create_sidebar
from components.viewer import create_viewer
from components.console import create_console

app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.CYBORG],
    suppress_callback_exceptions=True,
    title="CR3BP Explorer",
)

app.layout = dbc.Container(
    [
        # Shared state stores
        dcc.Store(id="system-params", data={}),
        dcc.Store(id="corrector-result", data=None),

        dbc.Row(
            [
                # Left sidebar
                create_sidebar(),

                # Main area: viewer + console
                dbc.Col(
                    [
                        create_viewer(),
                        create_console(),
                    ],
                    width=9,
                    className="p-0",
                ),
            ],
            className="g-0",
        ),
    ],
    fluid=True,
)

# Import callbacks to register them with the app
import callbacks.system   # noqa: F401
import callbacks.family   # noqa: F401
import callbacks.corrector  # noqa: F401

if __name__ == "__main__":
    app.run(debug=True, host="127.0.0.1", port=8050)
