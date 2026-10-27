# pctc in Lean

The test suite checks the numerics of `PctcLoop` and `Capacities` against closed forms and against
two independent SDP implementations. That is evidence, not proof: agreement to `1e-15` does not
show that an identity holds, and agreement with a solver over the whole Hermitian cone does not show
why a restricted search should match it. This directory proves the identities the Scala code relies
on, in Lean 4 with Mathlib.

## Scope

**The noiseless loop** (`Pctc/Loop.lean`). For a finite chronology-respecting system `S`, a finite
CTC system `K` of dimension `d`, any matrix `U` on `S ⊗ K` (it need not be unitary) and any matrix
`ρ` on `S`, the numerator of Eq. (1) of arXiv:2509.08965v3 for a noiseless loop (`N = id`), the map
`T(X) = U X Uᴴ` and a trivial reference system equals `C ρ Cᴴ / d²` with `C = Tr_K[U]`. After the
loop divides out its post-selection probability, the output is the Lloyd et al. evolution
`ρ ↦ C ρ Cᴴ / Tr[C ρ Cᴴ]` of arXiv:1005.2219. For demonstration (a), `U = RY(θ)` on the CTC qubit
alone, the post-selection probability is `cos²(θ/2)` and is exactly zero at `θ = π`.

**The block-diagonal restriction** (`Pctc/BlockDiagonal.lean`). With `λ` folded into the variable,
Eqs. (9) and (10) are `min Tr Y` over `Y ≥ 0` with `π ⊗ Y - J ≥ 0`, and `max Tr W` over Hermitian
`W` with `J - π ⊗ W ≥ 0`. For the erasure channel's Choi operator `J`, every feasible `Y` (or `W`)
has a block-diagonal feasible counterpart with the same trace: no coherence between the qubit block
and the flag `|e⟩`. So the block-diagonal qutrits `Bloch` sweeps reach every value the full programs
reach, for every `π` and every `p`. The argument is the one in the scaladoc of `Bloch.blockQutrit`,
with one simplification: averaging over `{I, D}`, `D = diag(1, 1, -1)`, is already enough, so the
whole phase group `I₂ + e^{iφ}|e⟩⟨e|` is not needed. The general form is `maxInformation_pinch` and
`doeblinInformation_pinch`, for any channel whose Choi operator is invariant under `I ⊗ D`; by
`choi_invariant` that holds whenever every Kraus operator maps into a single block.

Not covered: a noisy loop `N`, a non-trivial reference system `R`, a general bipartite `T`,
post-selected teleportation, Deutsch's and Allen's models, the closed forms for `I_max` and `I_doe`,
the grid and bisection that approximate them, and the restriction of the qubit sweep to a bounded
Bloch radius. Those remain checked numerically only.

## What to trust

A proof is only as good as its definitions. These are the ones to read against the Scala; the rest
is proof.

| Lean | Scala | Meaning |
| --- | --- | --- |
| `traceRight` | `DensityMatrix.partialTrace` | `Tr_K`, summing the right factor's diagonal |
| `phi` | `DensityMatrix.maximallyEntangled` | `Φ = \|Ω⟩⟨Ω\| / d`, `\|Ω⟩ = ∑ₖ \|k k⟩` |
| `loopOperator` | `PctcLoop.loopOperator` | `C = Tr_CTC[U]` |
| `loopInput` | `rhoRE.tensor(loopResource)` | `ρ ⊗ Φ_{BA'}` |
| `loopEvolved` | `t.applyTo(input, Seq(1, 2))` | `(U ⊗ I) (ρ ⊗ Φ) (U ⊗ I)ᴴ` |
| `loopWeighted` | `projector.multiply(evolved.matrix)` | `(I ⊗ Φ_{AA'}) · evolved` |
| `loopNumerator` | `weighted.partialTrace(Seq(2, 3))` | the numerator of Eq. (1) |
| `ry` | `StandardGate.RY` (Scotty) | `RY(θ)` |
| `erasureKraus` | `Channel.erasure` (Scotty) | the three Kraus operators, in Scotty's order |
| `erasureFlag` | basis element 2 | the flag `\|e⟩` |
| `erasureChoi` | `Capacities.choiState` | `∑ₖ (I ⊗ Kₖ) Φ (I ⊗ Kₖ)ᴴ` on `(A', B)` |
| the two constraint forms | `OperatorOrder.minLogLambda` | Eqs. (9) and (10) with `λ` folded in |

`erasureChoi_kept` and `erasureChoi_flag` pin the entries of `erasureChoi` to `(1 - p) Φ` on the
qubit block and `p π ⊗ |e⟩⟨e|` on the flag, so the invariance is not true of a degenerate
definition by accident.

## Statements

| Result | Lean declaration |
| --- | --- |
| numerator of Eq. (1) is `C ρ Cᴴ / d²` | `Pctc.loopNumerator_eq` |
| post-selection probability is `Tr[C ρ Cᴴ] / d²` | `Pctc.trace_loopNumerator` |
| normalised output is `C ρ Cᴴ / Tr[C ρ Cᴴ]` | `Pctc.loop_normalized_eq` |
| grandfather paradox: `P = cos²(θ/2)` | `Pctc.grandfather_postSelection` |
| grandfather paradox at `θ = π`: `P = 0` | `Pctc.grandfather_postSelection_pi` |
| `I_max`: pinching keeps feasibility and trace | `Pctc.maxInformation_pinch` |
| `I_doe`: pinching keeps feasibility and trace | `Pctc.doeblinInformation_pinch` |
| block-preserving Kraus operators give an invariant Choi operator | `Pctc.choi_invariant` |
| erasure: `I_max` optimum can be taken block diagonal | `Pctc.erasure_maxInformation_blockDiagonal` |
| erasure: `I_doe` optimum can be taken block diagonal | `Pctc.erasure_doeblinInformation_blockDiagonal` |

## Building

Lean `v4.35.0`, Mathlib `v4.35.0` (pinned in `lake-manifest.json`).

```bash
cd lean
lake exe cache get              # prebuilt Mathlib; without it Mathlib compiles from source
lake build
lake env lean check_axioms.lean
```

Every theorem should report `[propext, Classical.choice, Quot.sound]` and nothing else. `sorryAx`
anywhere would mean a gap.

## A second kernel

`kolla-check.sh` exports the dependency closure of every theorem in `check_axioms.lean` with
[lean4export](https://github.com/leanprover/lean4export) and re-checks it with
[kolla](https://github.com/oluies/kolla), a Lean kernel written independently in Scala, so the
verdict does not rest on Lean's own kernel alone.

```bash
LEAN4EXPORT=/path/to/lean4export KOLLA_JAR=/path/to/kolla.jar ./kolla-check.sh
```

At the time of writing the closure is 18,168 declarations (79 MB), all accepted with none carried
over, and kolla reports the same three axioms for each headline theorem. It takes about 15 seconds
on an M5 Pro after a 90-second export.
