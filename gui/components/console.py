"""
console.py

Bottom panel: terminal-style log output for corrector iterations.
"""

import dash_bootstrap_components as dbc
from dash import html, dcc


def create_console():
    return html.Div(
        [
            html.Div(
                [
                    html.Span("Console", className="console-title"),
                    dbc.Button(
                        "Clear", id="clear-console-btn", size="sm",
                        color="secondary", outline=True, className="console-clear-btn",
                    ),
                ],
                className="console-header",
            ),
            html.Pre(
                id="console-output",
                children="CR3BP Explorer ready.\n",
                className="console-pre",
            ),
            dcc.Interval(id="console-poll-interval", interval=200, disabled=True),
        ],
        className="console-panel",
    )
