"""
viewer.py

Center panel: 3D Plotly scatter plot for orbit visualization.
"""

import dash_bootstrap_components as dbc
from dash import html, dcc
import plotly.graph_objects as go


def _empty_figure():
    """Create an empty dark-themed 3D figure."""
    fig = go.Figure()
    fig.update_layout(
        template="plotly_dark",
        paper_bgcolor="#0a0a1a",
        plot_bgcolor="#0a0a1a",
        scene=dict(
            bgcolor="#0a0a1a",
            xaxis=dict(
                backgroundcolor="#0a0a1a",
                gridcolor="#1a1a3a",
                showbackground=True,
                title="X",
                color="#5070aa",
            ),
            yaxis=dict(
                backgroundcolor="#0a0a1a",
                gridcolor="#1a1a3a",
                showbackground=True,
                title="Y",
                color="#5070aa",
            ),
            zaxis=dict(
                backgroundcolor="#0a0a1a",
                gridcolor="#1a1a3a",
                showbackground=True,
                title="Z",
                color="#5070aa",
            ),
            aspectmode="data",
        ),
        margin=dict(l=0, r=0, t=30, b=0),
        font=dict(family="JetBrains Mono, monospace", color="#8090b0"),
        title=dict(text="CR3BP Explorer", font=dict(size=14, color="#7090ff")),
    )
    return fig


def create_viewer():
    return html.Div(
        [
            dcc.Loading(
                id="viewer-loading",
                type="circle",
                color="#7090ff",
                children=dcc.Graph(
                    id="orbit-3d-plot",
                    figure=_empty_figure(),
                    style={"height": "75vh"},
                    config={"scrollZoom": True},
                ),
            ),
        ],
        className="viewer-panel",
    )
