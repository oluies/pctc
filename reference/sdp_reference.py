# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy>=2", "cvxpy>=1.7"]
# ///
"""Independent reference values for the max- and Doeblin informations.

The Scala tests otherwise check themselves against closed forms transcribed from the same paper
they implement, so a transcription slip would pass on both sides. This solves the two semidefinite
programs of arXiv:2509.08965 Eqs. (9) and (10) directly from their definitions, with an off-the-shelf
interior-point solver, and writes the answers to a fixture the Scala suite reads back.

Both programs become linear once `lambda` is folded into the variable:

    I_max(N) = log2  min { tr Y : Y >= 0,          pi (x) Y >= N(Phi) }      with Y = lambda sigma
    I_doe(N) = -log2 max { tr W : W = W^dagger,    pi (x) W <= N(Phi) }      with W = tau / lambda

Nothing here is restricted to the block-diagonal or Bloch forms the Scala search sweeps, so this
also tests that those restrictions lose nothing.

Usage:  uv run reference/sdp_reference.py [--out PATH]
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import cvxpy as cp
import numpy as np

# The Doeblin program is unbounded exactly when its optimal trace is zero, and a solver reports
# that as a value a few times its own accuracy rather than as zero. The smallest legitimately finite
# trace anywhere in this sweep is 2^-8.57 = 0.0026, so anything under 1e-6 is unboundedness, not a
# small answer.
UNBOUNDED_BELOW = 1e-6


def depolarizing(p: float) -> list[np.ndarray]:
    """rho -> (1 - p) rho + p I/2."""
    x = np.array([[0, 1], [1, 0]], dtype=complex)
    y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    z = np.array([[1, 0], [0, -1]], dtype=complex)
    w = math.sqrt(p / 4)
    return [math.sqrt(1 - 3 * p / 4) * np.eye(2, dtype=complex), w * x, w * y, w * z]


def amplitude_damping(gamma: float) -> list[np.ndarray]:
    return [
        np.array([[1, 0], [0, math.sqrt(1 - gamma)]], dtype=complex),
        np.array([[0, math.sqrt(gamma)], [0, 0]], dtype=complex),
    ]


def erasure(p: float) -> list[np.ndarray]:
    """Qubit to qutrit; basis element 2 of the output is the erasure flag."""
    kept, lost = math.sqrt(1 - p), math.sqrt(p)
    return [
        kept * np.array([[1, 0], [0, 1], [0, 0]], dtype=complex),
        lost * np.array([[0, 0], [0, 0], [1, 0]], dtype=complex),
        lost * np.array([[0, 0], [0, 0], [0, 1]], dtype=complex),
    ]


CHANNELS = {
    "depolarizing": depolarizing,
    "amplitude damping": amplitude_damping,
    "erasure": erasure,
}


def choi(kraus: list[np.ndarray]) -> np.ndarray:
    """N(Phi_{A'A}) on A' (x) B, with Phi the normalised maximally entangled state."""
    d_in = kraus[0].shape[1]
    phi = np.zeros((d_in * d_in, 1), dtype=complex)
    for i in range(d_in):
        phi[i * d_in + i, 0] = 1 / math.sqrt(d_in)
    maximally_entangled = phi @ phi.conj().T

    identity = np.eye(d_in, dtype=complex)
    return sum(
        np.kron(identity, k) @ maximally_entangled @ np.kron(identity, k).conj().T for k in kraus
    )


def completeness_defect(kraus: list[np.ndarray]) -> float:
    total = sum(k.conj().T @ k for k in kraus)
    return float(np.max(np.abs(total - np.eye(total.shape[0]))))


def solve(objective: cp.Minimize | cp.Maximize, constraints: list) -> float | None:
    problem = cp.Problem(objective, constraints)
    for solver in (cp.CLARABEL, cp.SCS):
        try:
            problem.solve(solver=solver)
        except cp.error.SolverError:
            continue
        if problem.status in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
            return float(problem.value)
    return None


def max_information(j: np.ndarray, d_in: int) -> float:
    """log2 min { tr Y : Y >= 0, pi (x) Y >= J }."""
    out = j.shape[0] // d_in
    y = cp.Variable((out, out), hermitian=True)
    pi = np.eye(d_in, dtype=complex) / d_in
    value = solve(cp.Minimize(cp.real(cp.trace(y))), [y >> 0, cp.kron(pi, y) - j >> 0])
    return math.inf if value is None else math.log2(value)


def doeblin_information(j: np.ndarray, d_in: int) -> float:
    """-log2 max { tr W : W Hermitian, pi (x) W <= J }."""
    out = j.shape[0] // d_in
    w = cp.Variable((out, out), hermitian=True)
    pi = np.eye(d_in, dtype=complex) / d_in
    value = solve(cp.Maximize(cp.real(cp.trace(w))), [j - cp.kron(pi, w) >> 0])
    if value is None or value <= UNBOUNDED_BELOW:
        return math.inf
    return -math.log2(value)


def main() -> None:
    parser = argparse.ArgumentParser()
    resources = Path(__file__).parent.parent / "src/test/resources"
    parser.add_argument("--out", type=Path, default=resources)
    args = parser.parse_args()

    capacities = ["channel\tparameter\ti_max\ti_doe\tcompleteness_defect"]
    choi_entries = ["channel\tparameter\trow\tcol\tre\tim"]

    for name, build in CHANNELS.items():
        for step in range(11):
            parameter = step / 10
            kraus = build(parameter)
            j = choi(kraus)
            d_in = kraus[0].shape[1]

            capacities.append(
                f"{name}\t{parameter:.1f}\t{max_information(j, d_in):.12g}"
                f"\t{doeblin_information(j, d_in):.12g}\t{completeness_defect(kraus):.3e}"
            )

            for row in range(j.shape[0]):
                for col in range(j.shape[1]):
                    entry = j[row, col]
                    if abs(entry) > 1e-15:
                        choi_entries.append(
                            f"{name}\t{parameter:.1f}\t{row}\t{col}"
                            f"\t{entry.real:.17g}\t{entry.imag:.17g}"
                        )

    (args.out / "reference-capacities.tsv").write_text("\n".join(capacities) + "\n")
    (args.out / "reference-choi.tsv").write_text("\n".join(choi_entries) + "\n")

    print(f"wrote {len(capacities) - 1} capacity rows and {len(choi_entries) - 1} Choi entries to {args.out}")


if __name__ == "__main__":
    main()
