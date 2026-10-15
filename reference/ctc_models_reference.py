# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy>=2"]
# ///
"""Independent reference for both closed-timelike-curve models.

The Scala suite asserts the values published in J.-M. A. Allen, "Treating time travel quantum
mechanically", Phys. Rev. A 90, 042107 (2014), arXiv:1401.4933. This re-derives the same quantities
with numpy and a general eigensolver, so a mistake shared between the Scala code and my reading of
the paper would still show up as a mismatch.

Deutsch (arXiv:1401.4933 Eq. 1):  tau = Tr_CR[ U (rho_i (x) tau) U* ]
Post-selected (Eq. 9):            rho -> P rho P* / Tr[P rho P*],  P = Tr_CV(U)

Deutsch's condition is solved as an eigenproblem of the superoperator, never by iteration: for the
grandfather paradox the iteration has period two and never converges.

Usage:  uv run reference/ctc_models_reference.py
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

X = np.array([[0, 1], [1, 0]], dtype=complex)
H = np.array([[1, 1], [1, -1]], dtype=complex) / math.sqrt(2)


def bits_of(index: int, n: int) -> list[int]:
    return [(index >> (n - 1 - k)) & 1 for k in range(n)]


def index_of(bits: list[int]) -> int:
    return sum(b << (len(bits) - 1 - k) for k, b in enumerate(bits))


def controlled(n: int, control: int, target: int, gate: np.ndarray) -> np.ndarray:
    d = 2**n
    out = np.zeros((d, d), dtype=complex)
    for col in range(d):
        bits = bits_of(col, n)
        if bits[control] == 0:
            out[col, col] = 1
        else:
            for new in (0, 1):
                amp = gate[new, bits[target]]
                if amp != 0:
                    out[index_of(bits[:target] + [new] + bits[target + 1 :]), col] += amp
    return out


def on_qubit(n: int, target: int, gate: np.ndarray) -> np.ndarray:
    d = 2**n
    out = np.zeros((d, d), dtype=complex)
    for col in range(d):
        bits = bits_of(col, n)
        for new in (0, 1):
            amp = gate[new, bits[target]]
            if amp != 0:
                out[index_of(bits[:target] + [new] + bits[target + 1 :]), col] += amp
    return out


def swap(n: int, a: int, b: int) -> np.ndarray:
    d = 2**n
    out = np.zeros((d, d), dtype=complex)
    for col in range(d):
        bits = bits_of(col, n)
        bits[a], bits[b] = bits[b], bits[a]
        out[index_of(bits), col] = 1
    return out


def partial_trace(rho: np.ndarray, dims: list[int], keep: list[int]) -> np.ndarray:
    rho = rho.reshape(dims + dims)
    n = len(dims)
    for offset, axis in enumerate(sorted(set(range(n)) - set(keep), reverse=True)):
        rho = np.trace(rho, axis1=axis, axis2=axis + n - offset)
    kept = math.prod(dims[k] for k in keep)
    return rho.reshape(kept, kept)


def deutsch(u: np.ndarray, rho_i: np.ndarray, cr_dims: list[int], ctc_dim: int):
    """Returns (number of free parameters, maximum-entropy tau, induced output)."""
    dims = cr_dims + [ctc_dim]
    cr = math.prod(cr_dims)

    def step(tau: np.ndarray) -> np.ndarray:
        evolved = u @ np.kron(rho_i, tau) @ u.conj().T
        return partial_trace(evolved, dims, [len(cr_dims)])

    basis = [np.zeros((ctc_dim, ctc_dim), dtype=complex) for _ in range(ctc_dim**2)]
    for k, b in enumerate(basis):
        b[k // ctc_dim, k % ctc_dim] = 1
    super_op = np.column_stack([step(b).reshape(-1) for b in basis])

    vals, vecs = np.linalg.eig(super_op)
    fixed = [vecs[:, i].reshape(ctc_dim, ctc_dim) for i in range(len(vals)) if abs(vals[i] - 1) < 1e-8]

    hermitian = []
    for f in fixed:
        for part in (f + f.conj().T, -1j * (f - f.conj().T)):
            if np.max(np.abs(part)) > 1e-9:
                hermitian.append(part / 2)

    # Sweep the affine family on a grid and keep the positive member of greatest entropy.
    traced = next(h for h in hermitian if abs(np.trace(h)) > 1e-9)
    particular = traced / np.trace(traced)
    directions = [h - np.trace(h) * particular for h in hermitian if h is not traced]
    directions = [d / np.linalg.norm(d) for d in directions if np.linalg.norm(d) > 1e-9]

    # Independent directions only.
    independent: list[np.ndarray] = []
    for d in directions:
        r = d
        for b in independent:
            r = r - b * np.real(np.trace(b.conj().T @ r))
        if np.linalg.norm(r) > 1e-7:
            independent.append(r / np.linalg.norm(r))

    best_tau, best_entropy = particular, -1.0
    grid = np.linspace(-2, 2, 81)
    for coords in np.ndindex(*([len(grid)] * len(independent))):
        tau = particular + sum(grid[c] * independent[k] for k, c in enumerate(coords))
        w = np.linalg.eigvalsh(tau)
        if w.min() < -1e-9:
            continue
        s = float(-sum(x * math.log2(x) for x in w if x > 1e-12))
        if s > best_entropy:
            best_tau, best_entropy = tau, s

    evolved = u @ np.kron(rho_i, best_tau) @ u.conj().T
    return len(independent), best_tau, partial_trace(evolved, dims, list(range(len(cr_dims))))


def postselected(u: np.ndarray, rho_i: np.ndarray, cr: int, ctc_dim: int):
    """Returns (post-selection probability, normalised output or None)."""
    p = np.zeros((cr, cr), dtype=complex)
    for a in range(cr):
        for b in range(cr):
            for t in range(ctc_dim):
                p[a, b] += u[a * ctc_dim + t, b * ctc_dim + t]
    transformed = p @ rho_i @ p.conj().T
    weight = float(np.real(np.trace(transformed)))
    if weight <= 1e-12:
        return 0.0, None
    return weight / ctc_dim**2, transformed / weight


def purity(rho) -> float:
    return float(np.real(np.trace(rho @ rho)))


def entropy(rho) -> float:
    w = np.linalg.eigvalsh(rho)
    return float(-sum(x * math.log2(x) for x in w if x > 1e-12))


def main() -> None:
    ket = lambda *a: np.array(a, dtype=complex)
    dm = lambda v: np.outer(v, v.conj())
    r2 = 1 / math.sqrt(2)

    grandfather = controlled(2, 0, 1, X)
    unproven = swap(3, 1, 2) @ controlled(3, 0, 1, X) @ controlled(3, 2, 0, X)
    discriminator = controlled(2, 0, 1, H) @ swap(2, 0, 1)

    scenarios = [
        ("grandfather CNOT, traveller |1>", grandfather, dm(ket(0, 1)), [2], 2),
        ("grandfather CNOT, traveller |0>", grandfather, dm(ket(1, 0)), [2], 2),
        ("Ry(theta) on the curve, theta = 0.0 pi", on_qubit(2, 1, np.eye(2, dtype=complex)), dm(ket(1, 0)), [2], 2),
        ("Ry(theta) on the curve, theta = 0.5 pi",
         on_qubit(2, 1, np.array([[math.cos(math.pi/4), -math.sin(math.pi/4)],
                                  [math.sin(math.pi/4), math.cos(math.pi/4)]], dtype=complex)), dm(ket(1, 0)), [2], 2),
        ("Ry(theta) on the curve, theta = 1.0 pi",
         on_qubit(2, 1, np.array([[0.0, -1.0], [1.0, 0.0]], dtype=complex)), dm(ket(1, 0)), [2], 2),
        ("unproven theorem, |0>_B|0>_M", unproven, dm(ket(1, 0, 0, 0)), [2, 2], 2),
        ("discrimination, |0>", discriminator, dm(ket(1, 0)), [2], 2),
        ("discrimination, |->", discriminator, dm(ket(r2, -r2)), [2], 2),
        ("discrimination, |+>", discriminator, dm(ket(r2, r2)), [2], 2),
        ("discrimination, |1>", discriminator, dm(ket(0, 1)), [2], 2),
    ]

    rows = ["scenario\tdeutsch_parameters\tdeutsch_entropy\tdeutsch_output_purity\tpctc_probability\tpctc_output_purity"]
    for name, u, rho_i, cr_dims, ctc_dim in scenarios:
        params, tau, out = deutsch(u, rho_i, cr_dims, ctc_dim)
        prob, ps_out = postselected(u, rho_i, math.prod(cr_dims), ctc_dim)
        ps_purity = "nan" if ps_out is None else f"{purity(ps_out):.12g}"
        rows.append(
            f"{name}\t{params}\t{entropy(tau):.12g}\t{purity(out):.12g}\t{prob:.12g}\t{ps_purity}"
        )

    out_path = Path(__file__).parent.parent / "src/test/resources/reference-ctc-models.tsv"
    out_path.write_text("\n".join(rows) + "\n")
    print(f"wrote {len(rows) - 1} rows to {out_path}")
    print("\n".join(rows))


if __name__ == "__main__":
    main()
