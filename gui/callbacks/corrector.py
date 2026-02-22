"""
corrector.py

Callbacks for differential corrector execution and console output polling.
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from dash import Input, Output, State, callback, no_update, ctx
import plotly.graph_objects as go

from services.propagator import create_ode, propagate_orbit
from services.corrector_runner import CorrectorRunner
from astro import jacobi as compute_jacobi

# Module-level runner instance
_runner = CorrectorRunner()


def _parse_index_list(text):
    """Parse comma-separated index string into list of ints."""
    if not text or not text.strip():
        return []
    return [int(x.strip()) for x in text.split(",") if x.strip().isdigit()]


@callback(
    Output("console-poll-interval", "disabled", allow_duplicate=True),
    Output("console-output", "children", allow_duplicate=True),
    Output("corrector-result", "data", allow_duplicate=True),
    Input("correct-orbit-btn", "n_clicks"),
    State("ic-x", "value"),
    State("ic-y", "value"),
    State("ic-z", "value"),
    State("ic-vx", "value"),
    State("ic-vy", "value"),
    State("ic-vz", "value"),
    State("ic-period", "value"),
    State("fix-init-input", "value"),
    State("fix-end-input", "value"),
    State("system-params", "data"),
    State("console-output", "children"),
    prevent_initial_call=True,
)
def start_corrector(n_clicks, x, y, z, vx, vy, vz, period,
                    fix_init_str, fix_end_str, sys_params, console_text):
    if not sys_params:
        return True, (console_text or "") + "Error: No system selected.\n", no_update

    try:
        state_6 = [float(x), float(y), float(z), float(vx), float(vy), float(vz)]
        period = float(period)
    except (TypeError, ValueError):
        return True, (console_text or "") + "Error: Invalid input values.\n", no_update

    mu = sys_params["mu"]
    lstar = sys_params["lstar"]
    tstar = sys_params["tstar"]

    fix_init = _parse_index_list(fix_init_str)
    fix_end = _parse_index_list(fix_end_str)

    ode = create_ode(mu, lstar, tstar)

    msg = (
        f"--- Starting ASSET Correction ---\n"
        f"IC: [{x}, {y}, {z}, {vx}, {vy}, {vz}]\n"
        f"T = {period}, fix_init={fix_init}, fix_end={fix_end}\n"
    )

    _runner.start(ode, state_6, period, fix_init=fix_init, fix_end=fix_end)

    return False, (console_text or "") + msg, no_update


@callback(
    Output("console-output", "children", allow_duplicate=True),
    Output("console-poll-interval", "disabled", allow_duplicate=True),
    Output("corrector-result", "data", allow_duplicate=True),
    Input("console-poll-interval", "n_intervals"),
    State("console-output", "children"),
    State("system-params", "data"),
    prevent_initial_call=True,
)
def poll_corrector(n_intervals, console_text, sys_params):
    lines = _runner.poll_output()
    text = console_text or ""

    if lines:
        text += "\n".join(lines) + "\n"

    if _runner.is_done:
        if _runner.error:
            text += f"\nCorrection FAILED: {_runner.error}\n"
            return text, True, no_update

        traj = _runner.result
        cflag = _runner.convergence_flag

        if traj and len(traj) > 0:
            final_state = traj[0]
            mu = sys_params.get("mu", 0)
            jc = compute_jacobi(final_state[:6], mu)
            period_final = traj[-1][6] - traj[0][6]
            text += (
                f"\n--- Correction Complete ---\n"
                f"Convergence: {cflag}\n"
                f"State: [{final_state[0]:.10f}, {final_state[1]:.10f}, {final_state[2]:.10f}, "
                f"{final_state[3]:.10f}, {final_state[4]:.10f}, {final_state[5]:.10f}]\n"
                f"Period: {period_final:.10f}\n"
                f"Jacobi: {jc:.10f}\n"
            )
            # Store corrected trajectory for plotting
            result_data = {
                "traj": [[float(s[i]) for i in range(7)] for s in traj],
                "jacobi": float(jc),
                "period": float(period_final),
            }
            return text, True, result_data
        else:
            text += "\nCorrection returned empty trajectory.\n"
            return text, True, no_update

    return text, False, no_update


# --- Plot corrected orbit ---

@callback(
    Output("orbit-3d-plot", "figure", allow_duplicate=True),
    Input("corrector-result", "data"),
    State("orbit-3d-plot", "figure"),
    prevent_initial_call=True,
)
def plot_corrected_orbit(result_data, current_fig):
    if not result_data or "traj" not in result_data:
        return no_update

    traj = result_data["traj"]
    xs = [s[0] for s in traj]
    ys = [s[1] for s in traj]
    zs = [s[2] for s in traj]

    fig = go.Figure(current_fig)
    fig.add_trace(go.Scatter3d(
        x=xs, y=ys, z=zs,
        mode="lines",
        line=dict(color="#ffffff", width=4),
        name="Corrected",
        showlegend=True,
    ))
    # Starting point marker
    fig.add_trace(go.Scatter3d(
        x=[xs[0]], y=[ys[0]], z=[zs[0]],
        mode="markers",
        marker=dict(size=5, color="#00ff88", symbol="circle"),
        name="IC",
        showlegend=True,
    ))

    return fig


# --- Clear console ---

@callback(
    Output("console-output", "children"),
    Input("clear-console-btn", "n_clicks"),
    prevent_initial_call=True,
)
def clear_console(n_clicks):
    return "Console cleared.\n"
