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

**The closed forms** (`Pctc/ClosedForms.lean`). For the depolarizing, erasure and amplitude-damping
channels, built from Scotty's Kraus operators, every closed form in `ClosedForms.scala` for `I_max`
and `I_doe` is the optimal value of its program: `IsLeast (maxValues J) v` for Eq. (9) and
`IsGreatest (doeValues J) v` for Eq. (10). Each proof gives a feasible operator that reaches `v`
and a dual operator `Z = ∑ᵢ cᵢ vᵢ vᵢᴴ`, `cᵢ ≥ 0`, with `Tr_{A'}[Z] / d = I`, whose weak-duality
bound shows nothing feasible does better. The certificates are in the docstrings; for amplitude
damping, `I_max` is `(1 + k)² = 2 (1 + k) - γ` with `k = √(1 - γ)`, reached at
`Y = diag(1 + k, k² + k)`.

**The repetition code** (`Pctc/RepetitionCode.lean`), a side section outside the paper. For the
codeword `b, …, b` of length `2t + 1` and an error pattern `e`, the set of flipped positions,
majority vote recovers `b` whenever `#e ≤ t`, and some pattern of `t + 1` flips defeats it. With
every position flipped independently with probability `p`, the probability that it misreads the
codeword, summed over all patterns, is the binomial tail `∑_{m > n/2} C(n, m) pᵐ (1 - p)ⁿ⁻ᵐ` that
`RepetitionCode.phaseFlipFailureClosedForm` computes. The quantum side, that a phase flip on `|±⟩`
acts as a classical flip on the `X`-basis outcome, is checked numerically, not proved.

Not covered: a noisy loop `N`, a non-trivial reference system `R`, a general bipartite `T`,
post-selected teleportation, Deutsch's and Allen's models, the step from the folded programs back
to the `log2 λ` of Eqs. (9) and (10), the grid and bisection that approximate them, and the
restriction of the qubit sweep to a bounded Bloch radius. Those remain checked numerically only.

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
| `depolarizingKraus` | `Channel.depolarizing` (Scotty) | `√(1 - 3p/4) I`, `√(p/4) X`, `√(p/4) Y`, `√(p/4) Z` |
| `amplitudeDampingKraus` | `Channel.amplitudeDamping` (Scotty) | the two Kraus operators, in Scotty's order |
| `maxMixed` | `Capacities.maximallyMixed` | `π = I / d` |
| `maxValues`, `doeValues` | `OperatorOrder.minLogLambda` | the traces Eqs. (9) and (10) reach with `λ` folded in |

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
| weak duality for Eqs. (9) and (10) | `Pctc.trace_le_of_maxFeasible`, `Pctc.trace_ge_of_doeFeasible` |
| depolarizing: `I_max = log2(4(1 - p) + p)`, Eq. (S199) | `Pctc.depolarizing_maxInformation` |
| depolarizing: `I_doe = -log2 p`, Eq. (S198) | `Pctc.depolarizing_doeblinInformation` |
| erasure: `I_max = log2(4(1 - p) + p)`, Eq. (S216) | `Pctc.erasure_maxInformation` |
| erasure: `I_doe = -log2 p`, Eq. (S215) | `Pctc.erasure_doeblinInformation` |
| amplitude damping: `I_max = log2(2(1 + √(1 - γ)) - γ)`, Eq. (S223) | `Pctc.amplitudeDamping_maxInformation` |
| amplitude damping: `I_doe = -2 log2(1 - √(1 - γ))`, Eq. (S222) | `Pctc.amplitudeDamping_doeblinInformation` |
| majority vote on `2t + 1` copies corrects `t` flips | `Pctc.majority_corrects` |
| … and not `t + 1` | `Pctc.majority_fails` |
| its failure probability is the binomial tail | `Pctc.failureProbability_eq` |

`depolarizingChoi_eq`, `erasureChoi_eq` and `amplitudeDampingChoi_eq` pin each Choi operator
built from the Kraus operators to its explicit form, the same check `reference/sdp_reference.py`
makes entry by entry.

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

At the time of writing the closure of the 29 theorems in `check_axioms.lean` is 18,961
declarations, all accepted with none carried over, and kolla reports the same three axioms for
every one of them. The check itself takes seconds; the export and the `lake build` before it take
the time.
