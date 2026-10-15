# pctc — a post-selected closed timelike curve, simulated

A numerical model of the framework in
[arXiv:2509.08965](https://arxiv.org/abs/2509.08965), *Retrocausal capacity of a quantum channel:
Communicating through noisy closed timelike curves* (Kaiyuan Ji, Seth Lloyd, Mark M. Wilde), built on
the density-matrix layer of [Scotty](../scotty).

## What this program does, and what it does not

It **does** evaluate, with ordinary linear algebra on 2- to 24-dimensional complex matrices:

- the nonlinear map a noisy post-selected closed timelike curve (P-CTC) induces on a bipartite
  completely positive map — Eq. (1) of the paper;
- the special case of a unitary loop with a noiseless CTC, checked against
  `rho -> C rho C† / Tr[C rho C†]` with `C = Tr_CTC(U)`, the post-selected CTC evolution of
  [arXiv:1005.2219](https://arxiv.org/abs/1005.2219) (Lloyd, Maccone, Garcia-Patron, Giovannetti,
  Shikano);
- the grandfather paradox as a post-selection probability that vanishes;
- post-selected teleportation and the renormalisation the loop performs on it;
- the max-information `I_max`, the Doeblin information `I_doe`, and their sum `C_retro`, for the
  depolarizing, amplitude-damping and erasure channels, next to the ordinary classical capacity of
  the same channel.

It **does not** send information backward in time, and nothing in the output should be read as
evidence that anything can. There is no closed timelike curve here: there is a linear map, a
projection onto a maximally entangled state, and a division by the probability of that projection.
The physical claim the paper makes is conditional — *granted the existence of a noisy P-CTC*, these
are the information-theoretic limits — and this program only computes the right-hand side of that
conditional. It also does not run on quantum hardware, does not simulate more than a handful of
qubits, and does not use a semidefinite-programming solver (see
[The approximation](#the-approximation) below).

## Prerequisites

Scotty's `pctc` branch, installed locally:

```bash
cd ../scotty
git switch pctc
sbt publishLocal publishM2       # or: ./mill scotty.publishM2Local
```

This project is built and tested on the latest JDK (27 at the time of writing) and emits bytecode
for JDK 21.
It depends on `xyz.entangled::scotty:0.7.0-SNAPSHOT`; the only other non-test dependency is
`org.duckdb:duckdb_jdbc`, used only by the `Main` / `Results` / `Query` layer, with
`org.scalameta::munit` for tests.

There is no CI workflow here, because a hosted runner cannot resolve the locally published Scotty
snapshot without also checking out and publishing that repository.

## Running

```bash
sbt run                                   # all three demonstrations, written to results.duckdb
sbt "run --no-db"                         # print only
sbt "run --steps 13 --rounds 40"          # a finer search grid
sbt "runMain pctc.Query"                  # summarise what is stored

./mill pctc.run                           # the same via Mill
./mill pctc.runMain pctc.Query
```

Options:

| Option | Default | Meaning |
| --- | --- | --- |
| `--steps N` | 9 | points per axis, in the coarse grid and in each refinement scan |
| `--rounds N` | 30 | coordinate-refinement rounds; the search window halves each round |
| `--radius X` | 1000 | largest Bloch radius swept for unit-trace Hermitian operators |
| `--theta-steps N` | 20 | intervals in the grandfather sweep from 0 to pi, so `N+1` points |
| `--db PATH` | `results.duckdb` | database file to append to |
| `--no-db` | — | skip storage entirely |

## The three demonstrations

### (a) Grandfather paradox

The chronology-violating qubit is rotated by `RY(theta)` inside an otherwise empty loop. The loop
operator is `C = Tr[RY(theta)] = 2 cos(theta/2)`, so the post-selection probability is
`|Tr RY(theta) / d|² = cos²(theta/2)`:

```
   theta | theta/pi | P(post-select) | cos^2(theta/2) | abs error
---------+----------+----------------+----------------+----------
0.000000 |   0.0000 | 1.000000000000 | 1.000000000000 |  4.44e-16
0.785398 |   0.2500 | 0.853553390593 | 0.853553390593 |  4.44e-16
1.570796 |   0.5000 | 0.500000000000 | 0.500000000000 |  2.22e-16
2.356194 |   0.7500 | 0.146446609407 | 0.146446609407 |  8.33e-17
3.141593 |   1.0000 | 0.000000000000 | 0.000000000000 |  1.37e-48
```

At `theta = pi` the rotation is a bit flip: the qubit would have to come back as its own negation,
and the branch has probability exactly zero. The paradox is resolved by never happening.

### (b) Post-selected teleportation

A message qubit is Bell-measured against the future end of a maximally entangled pair and the `Phi+`
outcome is post-selected. That outcome has probability 1/4, and conditioned on it the past end
carries the message exactly, with no correction to apply. Running the same bipartite map inside the
loop of Eq. (1) produces the same 1/4 as the denominator, which the loop divides out:

```
stage                                       | P(post-select) |   renormalised | message fidelity
--------------------------------------------+----------------+----------------+-----------------
Bell measurement by hand, product message   | 0.250000000000 | 0.250000000000 |   1.000000000000
Bell measurement by hand, entangled message | 0.250000000000 | 0.250000000000 |   1.000000000000
P-CTC loop of Eq. (1), product message      | 0.250000000000 | 1.000000000000 |   1.000000000000
P-CTC loop of Eq. (1), entangled message    | 0.250000000000 | 1.000000000000 |   1.000000000000
```

Both a product message and half of an entangled pair are carried, with fidelity 1 in each case.

### (c) Retrocausal messaging through a noisy CTC

For each channel and each parameter in `{0, 0.1, …, 1}`:

```
I_max(N)  = inf { log2 lambda : N(Phi_{A'A}) <= lambda pi_{A'} (x) sigma_B }     Eq. (9)
I_doe(N)  = inf { log2 lambda : pi_{A'} (x) tau_B <= lambda N(Phi_{A'A}) }       Eq. (10)
C_retro   = I_max + I_doe                                                        Eq. (11)
```

with `sigma_B` a state and `tau_B` a unit-trace Hermitian operator that need not be positive. The
table prints `C_retro` next to the ordinary classical capacity of the same channel, using the
standard closed forms (cited in `ClosedForms.scala`): `1 - h₂(p/2)` for the depolarizing channel
(King 2003), `1 - p` for the erasure channel (Bennett–DiVincenzo–Smolin 1997), and the Holevo
information `max_t [h₂((1-gamma) t) - t h₂(gamma)]` for amplitude damping (Giovannetti–Fazio 2005).
A sample row, at `p = 0.5`:

```
channel           | param |    I_max |    I_doe |   C_retro | C_retro ref |    C_ord | beats C_ord | P(post-select)
------------------+-------+----------+----------+-----------+-------------+----------+-------------+---------------
depolarizing      |   0.5 | 1.321928 | 1.000000 |  2.321928 |    2.321928 | 0.188722 |         yes |       0.156250
```

`P(post-select)` is the largest post-selection probability a teleportation loop through this noisy
CTC can reach. By Lemma 1 of the paper the maximum achievable singlet fraction is `2^I_max / d_A²`,
and the loop contributes a further `1 / d_A²`, so the probability is `2^I_max / d_A⁴`. For a noiseless
CTC that is `4/16 = 1/4`, which is exactly what demonstration (b) measures.

`C_retro` is infinite at parameter 0, where the channel is the identity and the Doeblin information
diverges. That is the honest answer, not a numerical failure, and the paper's Figs. S8–S10 say the
same. At parameter 1 the channel carries nothing and `C_retro` ties the ordinary capacity at zero.
Everywhere in between `C_retro` exceeds it, which is the comparison the paper draws.

### Caveat on the regularisation

Eq. (11) is `C_retro(N) = I_max(N) + I_doe^inf(N)` with the *regularised* Doeblin information.
`I_max` is additive, so its regularisation is itself. `I_doe` is additive for the depolarizing
channel (mictodiactic covariance, Eqs. (S196)–(S198)) and for the erasure channel (Eq. (S215)), so
for those two the printed `C_retro` is the capacity. For amplitude damping `I_doe` is only a lower
bound on its regularisation, so there the printed `C_retro` is an achievable rate; the paper makes
the same point below its Eq. (S226).

## The approximation

Both information measures are semidefinite programs. Rather than add an SDP solver, this program
sweeps the single operand over a Bloch parametrisation and finds `lambda` by bisection for each
swept operator. Every value reported is therefore an **upper bound on the true infimum, approached
from above as the grid is refined**.

What is swept:

| Quantity | Operand | Parametrisation | Box |
| --- | --- | --- | --- |
| `I_max`, qubit output | state `sigma_B` | `(I + r n̂·σ)/2`, spherical | `r ∈ [0,1]`, `theta ∈ [0,pi]`, `phi ∈ [0,2pi]` |
| `I_max`, qutrit output | state `sigma_B` | block diagonal, qubit block of trace `w` plus flag `1-w` | `w ∈ [0,1]` plus the above |
| `I_doe`, qubit output | unit-trace Hermitian `tau_B` | `(I + r n̂·σ)/2` with `r = u/(1-u)` | `u ∈ [0, R/(1+R)]` plus the angles |
| `I_doe`, qutrit output | unit-trace Hermitian `tau_B` | block diagonal, block trace `t` | `t ∈ [-2,2]` plus the above |

Three things are worth knowing about this:

1. **The radius map matters.** The optimal `tau_B` for the amplitude-damping channel sits at Bloch
   radius `2(1 + sqrt(1-gamma))/gamma - 1`, which is 38 at `gamma = 0.1` and diverges as
   `gamma -> 0`. A uniform grid in Bloch coordinates cannot reach it; a uniform grid in
   `u = r/(1+r)` can. `--radius` sets the largest `r` the box covers.
2. **The qutrit restriction is not a loss.** For the erasure channel the swept operator is taken
   block diagonal. The channel's Choi operator commutes with `I ⊗ (I₂ + e^{i phi}|e><e|)` for every
   phase, both operator inequalities are linear in the swept operator, and averaging over that phase
   group preserves the trace as well as the inequality — so the optimum can be taken block diagonal.
3. **Infinity is detected structurally, within a stated range.** `log2 lambda` is bisected up to a
   ceiling of 32, with a floor of -40 below which the two operators are treated as equal; a pair
   that is infeasible at the ceiling is reported as infinite. The semidefinite
   test normalises by `lambda` rather than multiplying by it, so one relative tolerance
   (`1e-12`, against a Jacobi accuracy of about `1e-15`) stays meaningful across the whole range.
   Without that normalisation a genuinely infeasible inequality passes at a large `lambda` and
   infinity is silently reported as a finite 30-odd bits.

### Resolution actually achieved

`CapacitiesSuite` asserts agreement with the closed forms, and `Query` reports the convergence from
the stored rows. Largest absolute deviation over all finite reference values, in bits:

| steps | rounds | depolarizing `I_max` / `I_doe` | amplitude damping | erasure |
| --- | --- | --- | --- | --- |
| 5 | 4 | 5.4e-11 / 6.1e-11 | 8.2e-05 / 4.4e-02 | 2.3e-02 / 6.6e-11 |
| 7 | 10 | 5.4e-11 / 6.1e-11 | 6.9e-09 / 4.4e-04 | 2.0e-04 / 6.0e-11 |
| 7 | 20 | 5.4e-11 / 6.1e-11 | 6.8e-11 / 1.3e-07 | 2.3e-07 / 6.0e-11 |
| 9 | 30 | 5.4e-11 / 6.1e-11 | 6.7e-11 / 1.0e-09 | 2.1e-10 / 6.9e-11 |

The depolarizing optima are grid points (`sigma_B = tau_B = pi`), so they are exact at every
resolution; the floor of about `6e-11` is the bisection tolerance, not the grid.

### Cost

The sweep at `--steps 17 --rounds 50` visits 2,232,200 grid points and takes about 1.6 s. Three
things keep it there, none of which change a single digit of the output:

- The semidefinite test goes through `Cholesky.differenceIsPositiveDefinite` rather than computing
  eigenvalues. Only the sign of the smallest eigenvalue matters, a factorisation decides that
  without iterating, and fusing the difference into the factorisation avoids allocating two
  matrices of boxed complex numbers per test.
- A grid point that cannot beat the incumbent costs **one** semidefinite test. Only a point that
  improves on it pays for a bisection, and the coarse grid bisects to a loose 0.125 bits, since all
  it has to do is pick a cell. The refinement rounds run at full tolerance — they have to, because
  their window halves every round and a centre displaced by a loose comparison is one the later
  rounds can no longer walk back — and the winning point is re-evaluated once at full precision at
  the end.
- The coarse grid is scanned in parallel, each chunk carrying its own incumbent so the result stays
  deterministic.

For reference, the same sweep took 86 s before these three changes, and the worst deviation from
the closed forms is unchanged at `9.7e-10`.

## Results storage

Every computed row is appended to `results.duckdb`, one table per demonstration, stamped with a
`run_id`, a timestamp, the git commit of the code that produced it (suffixed `-dirty` for an unclean
tree), and the search resolution. Nothing is ever updated or deleted, so runs at different
resolutions sit side by side and can be compared.

```
grandfather    run_id, created_at, git_commit, demo, channel, parameter_name, parameter,
               post_selection_probability, expected_probability, absolute_error, search_resolution

teleportation  run_id, created_at, git_commit, demo, channel, stage,
               post_selection_probability, renormalized_probability, message_fidelity,
               search_resolution

retrocausal    run_id, created_at, git_commit, demo, channel, parameter_name, parameter,
               i_max, i_doe, c_retro,
               i_max_reference, i_doe_reference, c_retro_reference,
               ordinary_capacity, post_selection_probability,
               search_resolution, grid_steps, refinement_rounds, search_radius, evaluations
```

`Query` prints the two summaries worth having across runs — the parameter range over which
`C_retro` exceeds the ordinary capacity, and the convergence of `I_max` and `I_doe` as the grid is
refined:

```bash
sbt "runMain pctc.Query"
sbt "runMain pctc.Query --db /path/to/results.duckdb"
```

### Querying from the DuckDB CLI

```bash
duckdb results.duckdb
```

```sql
-- what is in here
SELECT run_id, min(created_at) AS at, any_value(search_resolution) AS grid, count(*) AS rows
FROM retrocausal GROUP BY run_id ORDER BY at;

-- where the retrocausal capacity beats the ordinary one, most recent run
SELECT channel, min(parameter) AS from_p, max(parameter) AS to_p, count(*) AS points
FROM retrocausal
WHERE run_id = (SELECT run_id FROM retrocausal ORDER BY created_at DESC LIMIT 1)
  AND (NOT isfinite(c_retro) OR c_retro > ordinary_capacity + 1e-6)
GROUP BY channel ORDER BY channel;

-- how the search converges as the grid is refined
SELECT channel, grid_steps, refinement_rounds,
       max(abs(i_max - i_max_reference)) FILTER (WHERE isfinite(i_max_reference)) AS i_max_error,
       max(abs(i_doe - i_doe_reference)) FILTER (WHERE isfinite(i_doe_reference)) AS i_doe_error
FROM retrocausal
GROUP BY ALL ORDER BY channel, grid_steps, refinement_rounds;

-- the full comparison table for one run
SELECT channel, parameter, i_max, i_doe, c_retro, ordinary_capacity, post_selection_probability
FROM retrocausal
WHERE run_id = (SELECT run_id FROM retrocausal ORDER BY created_at DESC LIMIT 1)
ORDER BY channel, parameter;

-- the grandfather paradox, worst deviation from cos^2(theta/2)
SELECT run_id, count(*) AS rows, max(absolute_error) AS worst FROM grandfather GROUP BY run_id;
```

## Layout

| File | Role |
| --- | --- |
| `PctcLoop.scala` | Eq. (1), and `Tr_CTC(U)` for the unitary special case |
| `Capacities.scala` | `I_max` and `I_doe` by search; the Choi state |
| `OperatorOrder.scala` | `X <= lambda Y` and the bisection for the smallest `lambda` |
| `Bloch.scala` | the swept operator families |
| `CoordinateSearch.scala` | coarse product grid plus coordinate refinement |
| `ClosedForms.scala` | reference values and ordinary classical capacities, with citations |
| `ChannelSpec.scala` | the three channel families and their references |
| `Demos.scala` | the three demonstrations, as pure computations |
| `Rows.scala` | the result row types |
| `TextTable.scala` | plain-text table rendering |
| `Main.scala` | printing and DuckDB writes |
| `Results.scala` | the DuckDB schema and inserts |
| `Query.scala` | the DuckDB summaries |

Everything above `Main`, `Results` and `Query` is free of any database dependency and prints nothing.

## References

- Kaiyuan Ji, Seth Lloyd, Mark M. Wilde, *Retrocausal capacity of a quantum channel: Communicating
  through noisy closed timelike curves*, [arXiv:2509.08965](https://arxiv.org/abs/2509.08965).
  Equation numbers in the source refer to v3; `Sxxx` numbers are from its Supplemental Material.
- Seth Lloyd, Lorenzo Maccone, Raul Garcia-Patron, Vittorio Giovannetti, Yutaka Shikano,
  *The quantum mechanics of time travel through post-selected teleportation*,
  [arXiv:1005.2219](https://arxiv.org/abs/1005.2219).
- Christopher King, *The capacity of the quantum depolarizing channel*, IEEE Trans. Inf. Theory 49,
  221 (2003).
- Charles H. Bennett, David P. DiVincenzo, John A. Smolin, *Capacities of quantum erasure channels*,
  Phys. Rev. Lett. 78, 3217 (1997).
- Vittorio Giovannetti, Rosario Fazio, *Information-capacity description of spin-chain correlations*,
  Phys. Rev. A 71, 032314 (2005).

## License

Apache 2.0, matching Scotty.
