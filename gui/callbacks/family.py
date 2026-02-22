"""
family.py

Callbacks for family browsing, orbit rendering, color-by, and click-to-populate.
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from dash import Input, Output, State, callback, html, no_update, ctx
import plotly.graph_objects as go

from services.data_loader import (
    get_available_families,
    get_available_modifiers,
    filter_orbits,
    orbit_row_to_ic,
)
from services.propagator import create_ode, propagate_orbit


# --- Cascade dropdowns ---

@callback(
    Output("family-dropdown", "options"),
    Output("family-dropdown", "value"),
    Input("lagrange-point-dropdown", "value"),
)
def update_family_options(lp):
    if not lp:
        return [], None
    families = get_available_families(lp)
    options = [{"label": f, "value": f} for f in families]
    return options, families[0] if families else None


@callback(
    Output("modifier-dropdown", "options"),
    Output("modifier-dropdown", "value"),
    Input("lagrange-point-dropdown", "value"),
    Input("family-dropdown", "value"),
)
def update_modifier_options(lp, family):
    if not lp or not family:
        return [], None
    mods = get_available_modifiers(lp, family)
    if mods == [""]:
        return [], None
    options = [{"label": m if m else "(none)", "value": m} for m in mods]
    default = mods[0] if mods else None
    return options, default


# --- Orbit data table ---

@callback(
    Output("orbit-data-table", "data"),
    Output("orbit-count-badge", "children"),
    Input("lagrange-point-dropdown", "value"),
    Input("family-dropdown", "value"),
    Input("modifier-dropdown", "value"),
)
def update_orbit_table(lp, family, modifier):
    if not lp or not family:
        return [], ""
    rows = filter_orbits(lp, family, modifier)
    table_data = []
    for r in rows:
        table_data.append({
            "x0": round(float(r["x0"]), 6),
            "z0": round(float(r["z0"]), 6),
            "yDot0": round(float(r["yDot0"]), 6),
            "T": round(float(r["T"]), 4),
            "JacobiValue": round(float(r["JacobiValue"]), 6),
            "Stability_Index": round(float(r["Stability_Index"]), 4),
        })
    badge = html.Span(f"{len(rows)} orbits", className="orbit-count")
    return table_data, badge


# --- Render family ---

@callback(
    Output("orbit-3d-plot", "figure"),
    Input("render-family-btn", "n_clicks"),
    State("lagrange-point-dropdown", "value"),
    State("family-dropdown", "value"),
    State("modifier-dropdown", "value"),
    State("system-params", "data"),
    State("color-by-dropdown", "value"),
    prevent_initial_call=True,
)
def render_family(n_clicks, lp, family, modifier, sys_params, color_by):
    if not lp or not family or not sys_params:
        return no_update

    mu = sys_params["mu"]
    lstar = sys_params["lstar"]
    tstar = sys_params["tstar"]
    lps_data = sys_params.get("lagrange_points", [])

    ode = create_ode(mu, lstar, tstar)
    rows = filter_orbits(lp, family, modifier)

    if not rows:
        return no_update

    # Propagate all orbits
    traces = []
    jacobis = []
    stab_indices = []
    periods = []

    for i, row in enumerate(rows):
        state_6, period, jacobi, si = orbit_row_to_ic(row)
        jacobis.append(jacobi)
        stab_indices.append(si)
        periods.append(period)

        traj = propagate_orbit(ode, state_6, period, n_points=300)
        xs = [s[0] for s in traj]
        ys = [s[1] for s in traj]
        zs = [s[2] for s in traj]

        traces.append((xs, ys, zs, i))

    # Color mapping
    if color_by == "Jacobi Constant":
        metric = jacobis
        colorbar_title = "Jacobi"
    elif color_by == "Stability Index":
        metric = stab_indices
        colorbar_title = "SI"
    elif color_by == "Period":
        metric = periods
        colorbar_title = "T"
    else:
        metric = None
        colorbar_title = None

    fig = go.Figure()

    if metric is not None:
        mn, mx = min(metric), max(metric)
        rng = mx - mn if mx != mn else 1.0
        import plotly.colors as pc
        plasma = pc.sequential.Plasma
        n_colors = len(plasma)

        for (xs, ys, zs, idx) in traces:
            norm = (metric[idx] - mn) / rng
            ci = int(norm * (n_colors - 1))
            ci = max(0, min(ci, n_colors - 1))
            color = plasma[ci]
            fig.add_trace(go.Scatter3d(
                x=xs, y=ys, z=zs,
                mode="lines",
                line=dict(color=color, width=2),
                name=f"Orbit {idx}",
                customdata=[idx],
                hovertemplate=f"{colorbar_title}={metric[idx]:.4f}<extra>Orbit {idx}</extra>",
                showlegend=False,
            ))

        # Invisible trace for colorbar
        fig.add_trace(go.Scatter3d(
            x=[None], y=[None], z=[None],
            mode="markers",
            marker=dict(
                size=0.001,
                color=[mn, mx],
                colorscale="Plasma",
                colorbar=dict(title=colorbar_title, thickness=15, len=0.5),
                showscale=True,
            ),
            showlegend=False,
            hoverinfo="skip",
        ))
    else:
        for (xs, ys, zs, idx) in traces:
            fig.add_trace(go.Scatter3d(
                x=xs, y=ys, z=zs,
                mode="lines",
                line=dict(color="#00ccff", width=2),
                name=f"Orbit {idx}",
                customdata=[idx],
                showlegend=False,
            ))

    # P1 and P2 markers
    fig.add_trace(go.Scatter3d(
        x=[-mu], y=[0], z=[0],
        mode="markers+text",
        marker=dict(size=6, color="#ffaa00"),
        text=["P1"], textposition="top center",
        textfont=dict(color="#ffaa00", size=10),
        showlegend=False, hoverinfo="text",
    ))
    fig.add_trace(go.Scatter3d(
        x=[1 - mu], y=[0], z=[0],
        mode="markers+text",
        marker=dict(size=4, color="#aaaaff"),
        text=["P2"], textposition="top center",
        textfont=dict(color="#aaaaff", size=10),
        showlegend=False, hoverinfo="text",
    ))

    # Lagrange point markers
    lp_names = ["L1", "L2", "L3", "L4", "L5"]
    if lps_data:
        for i, name in enumerate(lp_names):
            fig.add_trace(go.Scatter3d(
                x=[lps_data[i][0]], y=[lps_data[i][1]], z=[0],
                mode="markers+text",
                marker=dict(
                    size=4,
                    symbol="diamond",
                    color="#ff5555" if name == lp else "#555577",
                ),
                text=[name], textposition="top center",
                textfont=dict(
                    color="#ff5555" if name == lp else "#555577",
                    size=9,
                ),
                showlegend=False, hoverinfo="text",
            ))

    fig.update_layout(
        template="plotly_dark",
        paper_bgcolor="#0a0a1a",
        plot_bgcolor="#0a0a1a",
        scene=dict(
            bgcolor="#0a0a1a",
            xaxis=dict(backgroundcolor="#0a0a1a", gridcolor="#1a1a3a", showbackground=True, title="X", color="#5070aa"),
            yaxis=dict(backgroundcolor="#0a0a1a", gridcolor="#1a1a3a", showbackground=True, title="Y", color="#5070aa"),
            zaxis=dict(backgroundcolor="#0a0a1a", gridcolor="#1a1a3a", showbackground=True, title="Z", color="#5070aa"),
            aspectmode="data",
        ),
        margin=dict(l=0, r=0, t=30, b=0),
        font=dict(family="JetBrains Mono, monospace", color="#8090b0"),
        title=dict(text=f"{lp} {family}{' ' + modifier if modifier else ''} — {len(rows)} orbits", font=dict(size=14, color="#7090ff")),
    )

    return fig


# --- Click orbit table row to populate corrector inputs ---

@callback(
    Output("ic-x", "value"),
    Output("ic-y", "value"),
    Output("ic-z", "value"),
    Output("ic-vx", "value"),
    Output("ic-vy", "value"),
    Output("ic-vz", "value"),
    Output("ic-period", "value"),
    Input("orbit-data-table", "selected_rows"),
    State("orbit-data-table", "data"),
    State("lagrange-point-dropdown", "value"),
    State("family-dropdown", "value"),
    State("modifier-dropdown", "value"),
    prevent_initial_call=True,
)
def populate_corrector_from_table(selected_rows, table_data, lp, family, modifier):
    if not selected_rows or not table_data:
        return (no_update,) * 7

    row_idx = selected_rows[0]
    if row_idx >= len(table_data):
        return (no_update,) * 7

    # Get the full row from the original catalog
    rows = filter_orbits(lp, family, modifier)
    if row_idx >= len(rows):
        return (no_update,) * 7

    r = rows[row_idx]
    return (
        float(r["x0"]),
        float(r["y0"]),
        float(r["z0"]),
        float(r["xDot0"]),
        float(r["yDot0"]),
        float(r["zDot0"]),
        float(r["T"]),
    )
