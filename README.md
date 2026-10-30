# pctc — a post-selected closed timelike curve, simulated

**[Run it in your browser →](https://oluies.github.io/pctc/)** — the same computation, compiled to
JavaScript with Scala.js. No server does any of the arithmetic.

A numerical model of the framework in
[arXiv:2509.08965](https://arxiv.org/abs/2509.08965), *Retrocausal capacity of a quantum channel:
Communicating through noisy closed timelike curves* (Kaiyuan Ji, Seth Lloyd, Mark M. Wilde), built on
the density-matrix layer of [Scotty](https://github.com/oluies/scotty/releases/tag/v0.7.0) (release
0.7.0 of a fork of [entangled-xyz/scotty](https://github.com/entangled-xyz/scotty)).

The channel layer this depends on also has a hardware implementation in Chisel, verified against the
same Scotty code: [**brokkr**](https://oluies.github.io/brokkr/) — see
[below](#the-channel-layer-in-hardware).

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

Scotty 0.7.0, installed locally from its release tag:

```bash
git clone -b v0.7.0 https://github.com/oluies/scotty.git
cd scotty
sbt publishLocal publishM2       # or: ./mill scotty.publishM2Local
```

0.7.0 is also [published to GitHub Packages](https://github.com/oluies/scotty/releases/tag/v0.7.0);
resolving it from there needs a token with `read:packages`, as Scotty's README describes.

This project targets **Scala 3.9 LTS**, is built and tested on the latest JDK (27 at the time of
writing), and emits bytecode for JDK 21.
It depends on `io.github.oluies::scotty:0.7.0`; the only other non-test dependency is
`org.duckdb:duckdb_jdbc`, used only by the `Main` / `Results` / `Query` layer, with
`org.scalameta::munit` for tests.

CI builds Scotty from source before every job, since it is not on Maven Central and GitHub Packages
needs a token to read: `.github/workflows/ci.yml` checks out `oluies/scotty` at tag `v0.7.0`, runs `publishLocal publishM2`
there, and only then compiles and tests this project — under sbt on JDK 21, 25 and 27, and under
Mill on 27. A third job regenerates the cvxpy fixture and reports any drift without failing the
build.

## In the browser

`web/` compiles the computation to JavaScript and publishes it to GitHub Pages. It deliberately does
**not** restructure this project into a cross-build: the only platform-specific code is `Platform`,
which hands the coarse grid to a thread pool on the JVM and scans it in order in a browser, and
everything that cannot run in a browser — the DuckDB layer and the command-line entry point — is
excluded by filename in `build.sbt`. Scotty's sources are compiled from the sibling checkout, which
the JVM build already requires.

```bash
sbt site                       # assembles web/target/site (index.html + linked JavaScript)
sbt web/fastLinkJS             # quick link while iterating
```

Two things had to change to link for JavaScript, both improvements in their own right:

- `java.util.Locale` is not implemented in Scala.js at all, so `String.format(Locale.ROOT, ...)`
  could not link. `TextTable.number` now goes through `java.math.BigDecimal`, which is locale-free
  by construction and implemented on both platforms. It parses `value.toString` rather than the
  double, because `java.util.Formatter` rounds the *shortest decimal representation* and not the
  exact binary expansion — building `BigDecimal` from the double directly disagreed with the
  formatter on 21 of 5369 cases, and parsing `toString` agrees on all of them.
- `String.format`'s `%,d` grouping flag reaches `java.text.DecimalFormat`, also unimplemented;
  `TextTable.grouped` replaces it.

Scala.js cannot read classpath resources and has no DuckDB, so `ReferenceSuite` and `ResultsSuite`
cannot run there. `web/src/main/scala/pctc/SmokeMain.scala` checks the published closed forms of
arXiv:2509.08965 and the worked outputs of arXiv:1401.4933 against the linked JavaScript instead,
under `node`, and the Pages workflow fails the build on any mismatch. All 22 checks pass; the sweep
runs about 97,000 grid points in roughly 200 ms in a browser.

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

### (d) Three models of a closed timelike curve

The two standard quantum treatments plus Allen's, on the same circuits.

Deutsch's model (D. Deutsch, Phys. Rev. D **44**, 3197 (1991)) demands that the chronology-violating
system come back as it left:

```
tau = Tr_CR[ U (rho_i (x) tau) U* ]
```

`DeutschCtc` solves that as the kernel of `D - I`, never by iteration. Iterating does not generally
converge: for the grandfather paradox the map is `tau -> X tau X`, whose orbit has period two, so a
fixed-point loop oscillates forever and whatever it returns after a fixed number of steps is not a
solution at all. The kernel is exact, and its dimension reports Deutsch's uniqueness ambiguity
directly rather than hiding it; when the solution is not unique, Deutsch's evolutionary principle
picks the member of greatest von Neumann entropy.

Allen's own proposal, the transition-probability curve (T-CTC), averages over chronology-violating
pure states weighted by their transition probability — his Eq. (24). Carrying that integral out over
the Fubini-Study measure collapses it to Eq. (33), which is what `TctcLoop` implements:

```
rho_f = z [ C rho C* + Tr_CV( U (rho (x) 1) U* ) ],    C = Tr_CV(U)
```

The first term is exactly the post-selected numerator; the second is an ordinary quantum channel.
Allen reads it as "a weighted mixture of the corresponding P-CTC equation of motion with an ordinary
quantum channel".

The three diverge sharply. **Deutsch's condition is always satisfiable**, at the cost of mixing the
curve. **Post-selection gives a paradoxical history probability zero** — it simply never happens.
**Allen's never forbids anything either**, because the channel term cannot vanish, and it pays for
that in mixedness.

```
scenario                               | D sol.  | D purity | P(post-select) | P purity | T purity | T coh.
---------------------------------------+---------+----------+----------------+----------+----------+-------
grandfather CNOT, traveller |1>        | 1-param |   1.0000 |      forbidden |        - |   1.0000 | 0.0000
Ry(theta) on the curve, theta = 1.0 pi | 1-param |   1.0000 |      forbidden |        - |   1.0000 | 0.0000
unproven theorem, |0>_B|0>_M           | 1-param |   0.5000 |         0.5000 |   1.0000 |   0.6250 | 0.5000
discrimination, |0>                    | unique  |   1.0000 |         0.3750 |   1.0000 |   0.5918 | 0.4286
discrimination, |->                    | unique  |   1.0000 |         0.3750 |   1.0000 |   0.5918 | 0.4286
discrimination, |+>                    | unique  |   0.6667 |         0.1250 |   1.0000 |   0.5200 | 0.2000
```

`T coh.` is the share of the T-CTC output carried by its post-selected term. Two rows are worth
dwelling on. On the **grandfather** circuits `C|psi> = 0`, so that share is zero and the whole output
comes from the channel — which is exactly why post-selection forbids those histories and Allen's
model does not. On **discrimination**, `|0>` and `|->` give T-CTC *identical* outputs, so it cannot
tell them apart; Deutsch and post-selection both can, with certainty.

T-CTC purity is **not** bracketed by the other two. On the unproven theorem it sits between them
(0.625 against 0.5 and 1.0, which is Allen's Eq. (35)), but on discrimination it is more mixed than
either. An earlier version of this README claimed otherwise and a test caught it.

The circuits and their solutions are published, not invented here: the unproven-theorem and
state-discrimination circuits and their answers under *all three* models come from J.-M. A. Allen,
*Treating time travel quantum mechanically*, Phys. Rev. A **90**, 042107 (2014),
[arXiv:1401.4933](https://arxiv.org/abs/1401.4933). `DeutschCtcSuite` asserts those published values
directly — the one-parameter family `alpha|0><0| + (1-alpha)|1><1|`, the maximum-entropy choice
`alpha = 1/2`, the classically correlated Deutsch output against the pure entangled post-selected
one, the perfect discrimination of the non-orthogonal `|0>` and `|->`, and the loop operator
`P = |0><0| + |1><-|`, and Eq. (35)'s transition-probability output
`1/2|00><00| + 1/4|00><11| + 1/4|11><00| + 1/2|11><11|`, which `TctcLoopSuite` asserts to 1e-12.

The uncomfortable row is discrimination: a Deutsch curve tells two **non-orthogonal** states apart
with certainty. Allen built T-CTCs partly to avoid that, and the table shows it working.

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

### Choosing the grid

Measured, not guessed. Worst deviation from the closed forms over the whole sweep, against wall time
on an M5 Pro:

| steps | rounds | time | worst error |
| --- | --- | --- | --- |
| 9 | 20 | 0.22 s | 2.4e-07 |
| **9** | **30** | **0.30 s** | **9.7e-10** |
| 9 | 40 | 0.39 s | 9.7e-10 |
| 11 | 30 | 0.38 s | 9.7e-10 |
| 13 | 30 | 0.46 s | 9.7e-10 |
| 17 | 30 | 0.71 s | 9.7e-10 |
| 17 | 50 | 1.15 s | 9.7e-10 |

`9.7e-10` is the floor, set by the bisection tolerance rather than the grid. Three things follow.
The default of `steps=9, rounds=30` is the cheapest way to reach that floor, so raising either knob
buys nothing. Rounds past 30 are pure waste. And **rounds matter more than coarse resolution**:
`steps=25, rounds=15` is both slower *and* less accurate (2.9e-06) than `steps=17, rounds=30`.

Refinement rounds are not uniformly productive, which is why a "stop when a round stops improving"
rule is unsafe here: for the erasure channel the search improves at rounds 0 and 1, then nothing
until round 36. For the depolarizing channel the optimum is a coarse grid point and *every*
refinement round after the first is idle. `RoundsVsAccuracy` regenerates the table above.

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

Two things that were measured and found *not* to help, recorded so they are not retried blindly.
Thread count barely moves the sweep: the parallel coarse grid does scale (3.6x from one thread to
eighteen on its own), but refinement is about 77% of a search and coordinate descent is inherently
sequential across axes, so Amdahl caps the whole thing near 1.3x. And handing the semidefinite test
`pi (x) sigma` in block form, so the Kronecker product is never materialised, came out about 25%
*slower*: the product is built once per grid point, whereas the block indexing it replaces is paid
on every bisection step, several per point, and again for every entry of the magnitude scan.

The lever that did work was the representation. `ComplexMatrix` now stores interleaved primitives
rather than `Vector[Complex]`, and the kernels read them through `realAt`/`imagAt` instead of
allocating a `Complex` per entry. Measured back to back, with a hand-written primitive routine as an
unchanged control and discounting the 9% the control itself drifted:

| | boxed | primitive | control-adjusted |
| --- | --- | --- | --- |
| 6x6 positive-definiteness test | 0.185 us | 0.108 us | 1.57x |
| Kronecker product | 0.153 us | 0.056 us | 2.50x |
| `OperatorOrder.dominates` | 0.244 us | 0.138 us | 1.62x |
| sweep, steps=9 rounds=30 | 0.389 s | 0.307 s | **1.16x** |
| sweep, steps=17 rounds=50 | 1.207 s | 0.856 s | **1.29x** |

Accuracy and evaluation counts are identical either way. Note that the microbenchmark promised more
than it delivered: comparing the boxed path against hand-written primitive code suggested about 4x,
but `ComplexMatrix` cannot fully close that gap (bounds checks, call boundaries), and the
semidefinite test is not the whole of a grid point, so the sweep gained 1.2-1.3x rather than the 2x
extrapolated from the microbenchmark.

Three further directions were measured and rejected.

**Lazy precision.** The search's tolerances can be replaced by a bracket that refines only when a
comparison demands it, so precision is never configured per phase -- which is where this search once
had an accuracy bug. It works, and is a hair more accurate (9.722e-10 against 9.745e-10), but costs
about 1.5x: proving an ordering between two *equal* values forces refinement to the floor, and the
spherical Bloch parametrisation is degenerate at zero radius, so exact ties are common. Kept on the
`lazy-precision-search` branch rather than merged. It did expose something worth knowing: tie
breaking at that degenerate point is load bearing, and preferring the later of two tied points moves
the centre the refinement rounds start from, costing four orders of magnitude of accuracy.

**Cats Effect, or any other concurrency runtime.** There is nothing here for one to win. The
parallel coarse grid already scales 3.6x on its own, and refinement -- about 77% of a search -- is
coordinate descent, sequential across axes by construction. Finer-grained work distribution, the one
thing a work-stealing runtime would offer, was tested directly by lowering the chunk threshold from
512 points to 32: the sweep went from 0.185 s to 0.195 s, slightly worse. A fibre runtime would add
a dependency the project does not otherwise need, for a part of the work that is not the constraint.

**Scala Native.** The one real argument for it is warm-up: a single sweep is **47% JIT warm-up**
(0.349 s cold against 0.185 s steady state, measured by `WarmupBench`), and the command-line program
runs once and exits, so it pays that every time. Ahead-of-time compilation would recover most of it.
Against that, DuckDB's JDBC driver is a JVM library, so the entire `Results`/`Query` layer would have
to be replaced, and steady-state numeric throughput on tight double loops is where HotSpot is already
competitive with LLVM. Recovering 0.16 s per invocation does not pay for losing the storage layer.

A note on method, since this cost real time. Timings taken minutes apart on this machine drift by up
to 1.9x, enough to invent or erase any of the effects above. Every figure here comes from an A/B run
back to back in one sitting with an unchanged control measured alongside; two earlier "improvements"
in this file's history were machine drift and had to be retracted.

## Cross-validation against independent implementations

Most of the suite checks the numerics against closed forms transcribed from the same paper this
program implements, which would not catch a transcription slip. Two independent implementations in
`reference/` close that gap by building the three channels from their Kraus definitions and solving
Eqs. (9) and (10) as semidefinite programs, with no code, language, modelling layer or solver in
common with the Scala search or with each other:

| | stack | solver | accuracy on these programs |
| --- | --- | --- | --- |
| `reference/sdp_reference.py` | numpy + cvxpy | CLARABEL (interior point) | ~1e-8 |
| `reference/sdp_reference.jl` | Convex.jl | SCS (first order, tightened) | ~1e-9 |
| `reference/ctc_models_reference.py` | numpy | general eigensolver | ~1e-4 (grid) |

The third covers all three closed-timelike-curve models on the circuits of demonstration (d), so the
published values are checked twice: once as literal assertions against the paper, and once against
an implementation that shares no code with the Scala one.

Both programs become linear once `lambda` is folded into the variable:

```
I_max(N) = log2  min { tr Y : Y >= 0,     pi (x) Y >= N(Phi) }    with Y = lambda sigma
I_doe(N) = -log2 max { tr W : W = W',     pi (x) W <= N(Phi) }    with W = tau / lambda
```

This matters for more than arithmetic. The solvers optimise over the **whole Hermitian cone**, where
the Scala search sweeps a Bloch parametrisation and, for the erasure channel, only block-diagonal
operators. Agreement is evidence that neither restriction loses anything — the symmetry argument in
[The approximation](#the-approximation) is checked, not just asserted, and it is also
[proved](#proofs-in-lean). The Python fixture also pins
every Choi operator entry by entry, which tests the Kraus operators and the density-matrix layer
directly rather than through the capacities.

`ReferenceSuite` reads the committed fixtures from `src/test/resources`, so the suite needs neither
Python nor Julia. Regenerate them with:

```bash
uv run reference/sdp_reference.py
julia --project=reference reference/sdp_reference.jl
```

## Proofs in Lean

Cross-validation shows the numbers agree; it cannot show that an identity holds. `lean/` proves the
ones the unitary-loop code rests on, in Lean 4 with Mathlib: for a noiseless loop the numerator of
Eq. (1) is `C rho C† / d²` with `C = Tr_CTC(U)`, for any `U`; the normalised output is therefore the
Lloyd et al. evolution; and the grandfather paradox's post-selection probability is `cos²(theta/2)`,
exactly zero at `theta = pi`. The definitions follow `PctcLoop.evaluate` step by step.

It also proves the symmetry argument of [The approximation](#the-approximation): for the erasure
channel, every operator feasible in Eq. (9) or (10) has a block-diagonal feasible counterpart with
the same trace, so sweeping only block-diagonal qutrits loses nothing. That holds for every `p`, and
averaging over `{I, diag(1, 1, -1)}` is enough, so the full phase group isn't needed.

And it proves the six closed forms for `I_max` and `I_doe` the searches are checked against: for
each channel, built from Scotty's Kraus operators, the published value is the optimum of Eq. (9) or
(10), shown by a feasible operator that reaches it and a dual certificate that nothing beats it.
Outside the paper, it proves what the repetition-code side note rests on: majority vote on `2t + 1`
copies corrects every pattern of at most `t` flips and no more, and with independent flips it
fails with exactly the binomial tail probability the simulation is checked against. Recursive
majority on nine copies, a cheaper decoder, corrects only three flips where exact majority corrects
four.
[`lean/README.md`](lean/README.md) states the scope, the definitions a reader has to trust and what
is not covered (the noisy loop, the step back to `log2 lambda` and the grid search are not).

```bash
cd lean && lake exe cache get && lake build && lake env lean check_axioms.lean
```

`lean/kolla-check.sh` re-checks the proofs and everything they depend on in Mathlib with
[kolla](https://github.com/oluies/kolla), an independent Lean kernel, so the result doesn't rest on
Lean's kernel alone.

## The channel layer in hardware

[**brokkr**](https://oluies.github.io/brokkr/)
([oluies/brokkr](https://github.com/oluies/brokkr)) is a design assessment of whether a CUTE-style
Chisel tensor engine can accelerate quantum simulation, and its conclusion lands on the layer this
project is built from. A state vector offers too little arithmetic per byte to keep a wide multiplier
array busy; a density matrix offers 683 flop/byte at ten qubits against 15.9 for the best fused state
vector, because a channel is two complex matrix products per Kraus operator.

So the engine it builds, `ChannelEngine`, computes exactly `rho -> sum_k K_k rho K_k^dagger` —
`Channel.apply`, in gates. It is checked against Scotty's own result on 24 cases: four channels
(depolarizing, amplitude damping, dephasing and **erasure**) over several random states each, square
and rectangular, with trace preservation and Hermiticity asserted as well, since those hold whatever
the reference says.

Erasure matters here specifically. It maps a qubit to a qutrit, so its Kraus operators are `3 x 2`
and `rho`, `K` and the result are three different shapes — and it is one of the three channels whose
retrocausal capacity this program computes. The engine carried one square dimension for a while and
could not express it at all; that was recorded as a gap on this project's workload until it was
closed.

**This is not a cross-validation of the numbers above**, and should not be read as one. The engine
works in 24-bit fixed point (Q2.22) at dimensions 4 to 12, where its worst deviation from Scotty is
`1.1e-06`; the reference implementations in the previous section agree with the Scala search to
`1e-08` or better at up to 24 dimensions in double precision. It is an independent implementation of
the channel application, three levels down, not an independent check of the capacities.

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

comparison     run_id, created_at, git_commit, demo, scenario, source,
               deutsch_parameters, deutsch_entropy, deutsch_output_purity,
               pctc_probability, pctc_output_purity, pctc_forbids,
               tctc_output_purity, tctc_coherent_fraction, search_resolution

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
| `DeutschCtc.scala` | Deutsch's self-consistency condition, solved as a kernel |
| `TctcLoop.scala` | Allen's transition-probability curve, Eq. (33) |
| `CtcComparison.scala` | the two models side by side on published circuits |
| `Circuits.scala` | full matrices for multi-qubit gates |
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
- John-Mark A. Allen, *Treating time travel quantum mechanically*, Phys. Rev. A **90**, 042107
  (2014), [arXiv:1401.4933](https://arxiv.org/abs/1401.4933). Source of the D-CTC, P-CTC and T-CTC
  comparison circuits and their published solutions.
- David Deutsch, *Quantum mechanics near closed timelike lines*, Phys. Rev. D **44**, 3197 (1991).
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
