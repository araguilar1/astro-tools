"""
corrector_runner.py

Background-threaded differential corrector using ASSET's solvePeriodic.

Classes
-------
CorrectorRunner : Threaded corrector with stdout capture and queue-based output
"""

import sys
import os
import io
import threading
import queue
import traceback

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from CR3BP_Utils.solve_periodic_orbit import solvePeriodic


class CorrectorRunner:
    """
    Runs solvePeriodic in a background thread, capturing PSIOPT output.

    Usage
    -----
        runner = CorrectorRunner()
        runner.start(ode, state_6, period)
        while not runner.is_done:
            lines = runner.poll_output()
            ...
        result = runner.result  # corrected trajectory
    """

    def __init__(self):
        self._thread = None
        self._queue = queue.Queue()
        self._done = False
        self._result = None
        self._error = None
        self._convergence_flag = None

    @property
    def is_done(self):
        return self._done

    @property
    def result(self):
        return self._result

    @property
    def error(self):
        return self._error

    @property
    def convergence_flag(self):
        return self._convergence_flag

    def poll_output(self):
        """Drain all available lines from the output queue."""
        lines = []
        while True:
            try:
                line = self._queue.get_nowait()
                lines.append(line)
            except queue.Empty:
                break
        return lines

    def start(self, ode, state_6, period, fix_init=None, fix_end=None):
        """
        Launch the corrector in a background thread.

        Parameters
        ----------
        ode : CR3BP
            ASSET ODE instance
        state_6 : list of float
            [x, y, z, vx, vy, vz]
        period : float
            Full period guess (non-dimensional)
        fix_init : list of int, optional
            Indices to fix at front (default: [0,1,2,3,6])
        fix_end : list of int, optional
            Indices to fix at back (default: [1,3,5])
        """
        if fix_init is None:
            fix_init = [0, 1, 2, 3, 6]
        if fix_end is None:
            fix_end = [1, 3, 5]

        self._done = False
        self._result = None
        self._error = None
        self._convergence_flag = None

        def _run():
            captured = io.StringIO()
            old_stdout = sys.stdout
            try:
                sys.stdout = captured
                ig = list(state_6) + [0.0]
                traj, cflag = solvePeriodic(
                    ig, period, ode,
                    print_level=0,
                    fix_init=fix_init,
                    fix_end=fix_end,
                )
                sys.stdout = old_stdout

                output_text = captured.getvalue()
                for line in output_text.splitlines():
                    self._queue.put(line)

                self._result = traj
                self._convergence_flag = cflag
            except Exception as e:
                sys.stdout = old_stdout
                output_text = captured.getvalue()
                for line in output_text.splitlines():
                    self._queue.put(line)
                self._error = f"{type(e).__name__}: {e}\n{traceback.format_exc()}"
                self._queue.put(f"ERROR: {self._error}")
            finally:
                sys.stdout = old_stdout
                self._done = True

        self._thread = threading.Thread(target=_run, daemon=True)
        self._thread.start()
