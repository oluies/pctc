import Pctc.BlockDiagonal

/-!
# The closed forms for `I_max` and `I_doe`

`ClosedForms.scala` transcribes six values from the Supplemental Material of arXiv:2509.08965v3,
and the test suite checks the searches against them. This file proves that each one is the optimal
value of its semidefinite program.

With `λ` folded into the variable, as in `Pctc/BlockDiagonal.lean`, the programs are

* Eq. (9), `I_max = log2 v`: `v` is the least `Tr Y` over `Y ≥ 0` with `π ⊗ Y - J ≥ 0`;
* Eq. (10), `I_doe = -log2 v`: `v` is the greatest `Tr W` over Hermitian `W` with `J - π ⊗ W ≥ 0`;

with `J` the Choi operator on `(A', B)` and `π = I / d` maximally mixed on `A'`. `maxValues J` and
`doeValues J` are the sets of traces the feasible operators reach, and each result is an
`IsLeast` or `IsGreatest` on one of them: the value is attained, and no feasible operator does
better.

Each proof has two certificates. A feasible operator with the right trace shows the value is
reached. A dual operator `Z = ∑ᵢ cᵢ vᵢ vᵢᴴ`, `cᵢ ≥ 0`, with `Tr_{A'}[Z] / d = I` shows it cannot be
beaten: for feasible `Y`, `0 ≤ Tr[Z (π ⊗ Y - J)] = Tr Y - Tr[Z J]` (`trace_le_of_maxFeasible`), and
for feasible `W`, `0 ≤ Tr[Z (J - π ⊗ W)] = Tr[Z J] - Tr W` (`trace_ge_of_doeFeasible`).

The step from these folded programs to the `log2 λ` of Eqs. (9) and (10) is not formalised here.
-/

namespace Pctc

open Matrix
open scoped Kronecker ComplexOrder

noncomputable section

/-! ## The two programs -/

variable {A B : Type*} [Fintype A] [DecidableEq A] [Fintype B] [DecidableEq B]

/-- `π = I / d`, `Capacities.maximallyMixed`. -/
def maxMixed : Matrix A A ℂ := (Fintype.card A : ℂ)⁻¹ • 1

/-- The traces reached by operators feasible in Eq. (9), with `λ` folded into `Y = λ σ`. -/
def maxValues (J : Matrix (A × B) (A × B) ℂ) : Set ℝ :=
  {t | ∃ Y : Matrix B B ℂ, Y.PosSemidef ∧ (maxMixed ⊗ₖ Y - J).PosSemidef ∧ Y.trace = t}

/-- The traces reached by operators feasible in Eq. (10), with `λ` folded into `W = τ / λ`. -/
def doeValues (J : Matrix (A × B) (A × B) ℂ) : Set ℝ :=
  {t | ∃ W : Matrix B B ℂ, W.IsHermitian ∧ (J - maxMixed ⊗ₖ W).PosSemidef ∧ W.trace = t}

/-! ## Weak duality -/

/-- The partial trace over the left tensor factor, `Tr_{A'}`. -/
def traceLeft (Z : Matrix (A × B) (A × B) ℂ) : Matrix B B ℂ := fun b b' => ∑ a, Z (a, b) (a, b')

/-- `∑ᵢ cᵢ vᵢ vᵢᴴ`, positive semidefinite whenever every `cᵢ ≥ 0`. -/
def rankOneSum {n ι : Type*} [Fintype ι] (c : ι → ℝ) (v : ι → n → ℂ) : Matrix n n ℂ :=
  ∑ i, (c i : ℂ) • vecMulVec (v i) (star (v i))

theorem posSemidef_rankOneSum {n ι : Type*} [Fintype n] [Fintype ι] {c : ι → ℝ} (v : ι → n → ℂ)
    (hc : ∀ i, 0 ≤ c i) : (rankOneSum c v).PosSemidef :=
  posSemidef_sum _ fun i _ =>
    (posSemidef_vecMulVec_self_star (v i)).smul (Complex.zero_le_real.mpr (hc i))

theorem isHermitian_rankOneSum {n ι : Type*} [Fintype ι] (c : ι → ℝ) (v : ι → n → ℂ) :
    (rankOneSum c v).IsHermitian := by
  ext i j
  simp only [conjTranspose_apply, rankOneSum, Matrix.sum_apply, Matrix.smul_apply, vecMulVec_apply,
    Pi.star_apply, star_sum, smul_eq_mul, star_mul', Complex.star_def, Complex.conj_ofReal,
    Complex.conj_conj]
  refine Finset.sum_congr rfl fun _ _ => ?_
  ring

/-- A diagonal matrix with real entries. -/
def realDiagonal {n : Type*} [DecidableEq n] (d : n → ℝ) : Matrix n n ℂ := diagonal fun b => (d b : ℂ)

theorem posSemidef_realDiagonal {n : Type*} [Fintype n] [DecidableEq n] {d : n → ℝ}
    (h : ∀ b, 0 ≤ d b) : (realDiagonal d).PosSemidef :=
  PosSemidef.diagonal fun b => Complex.zero_le_real.mpr (h b)

theorem isHermitian_realDiagonal {n : Type*} [DecidableEq n] (d : n → ℝ) :
    (realDiagonal d).IsHermitian := by
  rw [IsHermitian, realDiagonal, diagonal_conjTranspose]
  congr 1
  ext b
  simp [Complex.conj_ofReal]

theorem trace_vecMulVec_mul {n : Type*} [Fintype n] (u : n → ℂ) (M : Matrix n n ℂ) :
    (vecMulVec u (star u) * M).trace = star u ⬝ᵥ (M *ᵥ u) := by
  simp only [Matrix.trace, diag_apply, mul_apply, vecMulVec_apply, dotProduct, mulVec,
    Finset.mul_sum, Pi.star_apply]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun _ _ => Finset.sum_congr rfl fun _ _ => ?_
  ring

theorem trace_rankOneSum_mul_nonneg {n ι : Type*} [Fintype n] [Fintype ι] {c : ι → ℝ}
    (v : ι → n → ℂ) (hc : ∀ i, 0 ≤ c i) {M : Matrix n n ℂ} (hM : M.PosSemidef) :
    0 ≤ (rankOneSum c v * M).trace := by
  rw [rankOneSum, Finset.sum_mul, trace_sum]
  refine Finset.sum_nonneg fun i _ => ?_
  rw [smul_mul_assoc, trace_smul, trace_vecMulVec_mul, smul_eq_mul]
  exact mul_nonneg (Complex.zero_le_real.mpr (hc i)) (hM.dotProduct_mulVec_nonneg _)

theorem trace_mul_maxMixed_kronecker (Z : Matrix (A × B) (A × B) ℂ) (Y : Matrix B B ℂ) :
    (Z * (maxMixed ⊗ₖ Y)).trace = (((Fintype.card A : ℂ)⁻¹ • traceLeft Z) * Y).trace := by
  simp [Matrix.trace, mul_apply, maxMixed, kroneckerMap_apply, one_apply, Fintype.sum_prod_type,
    traceLeft, Finset.mul_sum, Finset.sum_mul, ite_mul, mul_ite]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun _ _ => ?_
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun _ _ => Finset.sum_congr rfl fun _ _ => ?_
  ring

/-- Weak duality for Eq. (9): a dual certificate bounds every feasible trace from below. -/
theorem trace_le_of_maxFeasible {ι : Type*} [Fintype ι] {c : ι → ℝ} (v : ι → A × B → ℂ)
    (hc : ∀ i, 0 ≤ c i) (hZ : (Fintype.card A : ℂ)⁻¹ • traceLeft (rankOneSum c v) = 1)
    {J : Matrix (A × B) (A × B) ℂ} {Y : Matrix B B ℂ} (h : (maxMixed ⊗ₖ Y - J).PosSemidef) :
    (rankOneSum c v * J).trace ≤ Y.trace := by
  have := trace_rankOneSum_mul_nonneg v hc h
  rwa [mul_sub, trace_sub, trace_mul_maxMixed_kronecker, hZ, one_mul, sub_nonneg] at this

/-- Weak duality for Eq. (10): a dual certificate bounds every feasible trace from above. -/
theorem trace_ge_of_doeFeasible {ι : Type*} [Fintype ι] {c : ι → ℝ} (v : ι → A × B → ℂ)
    (hc : ∀ i, 0 ≤ c i) (hZ : (Fintype.card A : ℂ)⁻¹ • traceLeft (rankOneSum c v) = 1)
    {J : Matrix (A × B) (A × B) ℂ} {W : Matrix B B ℂ} (h : (J - maxMixed ⊗ₖ W).PosSemidef) :
    W.trace ≤ (rankOneSum c v * J).trace := by
  have := trace_rankOneSum_mul_nonneg v hc h
  rwa [mul_sub, trace_sub, trace_mul_maxMixed_kronecker, hZ, one_mul, sub_nonneg] at this

/-- A feasible operator of trace `t` and a dual certificate of value `t` make `t` the optimum of
Eq. (9). -/
theorem isLeast_maxValues {J : Matrix (A × B) (A × B) ℂ} {t : ℝ} (Y : Matrix B B ℂ)
    (hY : Y.PosSemidef) (hfeas : (maxMixed ⊗ₖ Y - J).PosSemidef) (htr : Y.trace = t)
    {ι : Type*} [Fintype ι] {c : ι → ℝ} (v : ι → A × B → ℂ) (hc : ∀ i, 0 ≤ c i)
    (hZ : (Fintype.card A : ℂ)⁻¹ • traceLeft (rankOneSum c v) = 1)
    (hZJ : (rankOneSum c v * J).trace = t) : IsLeast (maxValues J) t := by
  refine ⟨⟨Y, hY, hfeas, htr⟩, fun s ⟨Y', _, h', hs⟩ => ?_⟩
  have := trace_le_of_maxFeasible v hc hZ h'
  rw [hZJ, hs] at this
  exact_mod_cast this

/-- A feasible operator of trace `t` and a dual certificate of value `t` make `t` the optimum of
Eq. (10). -/
theorem isGreatest_doeValues {J : Matrix (A × B) (A × B) ℂ} {t : ℝ} (W : Matrix B B ℂ)
    (hW : W.IsHermitian) (hfeas : (J - maxMixed ⊗ₖ W).PosSemidef) (htr : W.trace = t)
    {ι : Type*} [Fintype ι] {c : ι → ℝ} (v : ι → A × B → ℂ) (hc : ∀ i, 0 ≤ c i)
    (hZ : (Fintype.card A : ℂ)⁻¹ • traceLeft (rankOneSum c v) = 1)
    (hZJ : (rankOneSum c v * J).trace = t) : IsGreatest (doeValues J) t := by
  refine ⟨⟨W, hW, hfeas, htr⟩, fun s ⟨W', _, h', hs⟩ => ?_⟩
  have := trace_ge_of_doeFeasible v hc hZ h'
  rw [hZJ, hs] at this
  exact_mod_cast this

/-! ## Vectors on `A' ⊗ B` for a qubit `A'` -/

variable {B' : Type*}

/-- The vector on `(A', B)` whose amplitude on `|a b⟩` is `m a b`. -/
def vec (m : Matrix (Fin 2) B' ℂ) : Fin 2 × B' → ℂ := fun p => m p.1 p.2

/-- `|Ω⟩ = |00⟩ + |11⟩`, so `Φ = |Ω⟩⟨Ω| / 2`. -/
def omega : Fin 2 × Fin 2 → ℂ := vec !![1, 0; 0, 1]
/-- `|00⟩ - |11⟩`, `|01⟩ + |10⟩`, `|01⟩ - |10⟩`: the other three Bell vectors, unnormalised. -/
def phiMinus : Fin 2 × Fin 2 → ℂ := vec !![1, 0; 0, -1]
def psiPlus : Fin 2 × Fin 2 → ℂ := vec !![0, 1; 1, 0]
def psiMinus : Fin 2 × Fin 2 → ℂ := vec !![0, 1; -1, 0]

/-- The same four vectors inside the qubit block of a qutrit `B`, and `|0 e⟩`, `|1 e⟩`. -/
def omegaQ : Fin 2 × Fin 3 → ℂ := vec !![1, 0, 0; 0, 1, 0]
def phiMinusQ : Fin 2 × Fin 3 → ℂ := vec !![1, 0, 0; 0, -1, 0]
def psiPlusQ : Fin 2 × Fin 3 → ℂ := vec !![0, 1, 0; 1, 0, 0]
def psiMinusQ : Fin 2 × Fin 3 → ℂ := vec !![0, 1, 0; -1, 0, 0]
def flag0 : Fin 2 × Fin 3 → ℂ := vec !![0, 0, 1; 0, 0, 0]
def flag1 : Fin 2 × Fin 3 → ℂ := vec !![0, 0, 0; 0, 0, 1]

/-- `|a b⟩`. -/
def ket (a b : Fin 2) : Fin 2 × Fin 2 → ℂ := fun p => if p = (a, b) then 1 else 0

theorem isHermitian_real_smul_one {n : Type*} [DecidableEq n] (r : ℝ) :
    ((r : ℂ) • (1 : Matrix n n ℂ)).IsHermitian := by
  rw [IsHermitian, conjTranspose_smul, conjTranspose_one, Complex.star_def, Complex.conj_ofReal]

/-- Expands certificates on `(A', B)` to entries. -/
macro "pctc_simp" : tactic =>
  `(tactic| simp [rankOneSum, vecMulVec_apply, vec, omega, phiMinus, psiPlus, psiMinus, ket,
    omegaQ, phiMinusQ, psiPlusQ, psiMinusQ, flag0, flag1,
    Fin.sum_univ_succ, Fintype.sum_prod_type, mul_apply, Matrix.trace, traceLeft, maxMixed,
    kroneckerMap_apply, one_apply, Matrix.smul_apply, realDiagonal, diagonal_apply])

/-- Closes what `pctc_simp` leaves: a ring identity, or arithmetic on literals. -/
macro "pctc_finish" : tactic => `(tactic| ((try ring_nf) <;> (try norm_num)))

/-! ## The depolarizing channel -/

/-- The Kraus operators of `Channel.depolarizing(p)`, in Scotty's order: `I`, `X`, `Y`, `Z`. -/
def depolarizingKraus (p : ℝ) : Fin 4 → Matrix (Fin 2) (Fin 2) ℂ
  | 0 => (Real.sqrt (1 - 3 * p / 4) : ℂ) • !![1, 0; 0, 1]
  | 1 => (Real.sqrt (p / 4) : ℂ) • !![0, 1; 1, 0]
  | 2 => (Real.sqrt (p / 4) : ℂ) • !![0, -Complex.I; Complex.I, 0]
  | 3 => (Real.sqrt (p / 4) : ℂ) • !![1, 0; 0, -1]

/-- `Capacities.choiState(Channel.depolarizing(p))`. -/
def depolarizingChoi (p : ℝ) : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ :=
  ∑ k, krausTerm (depolarizingKraus p k)

/-- `J = (1 - p) Φ + p I / 4`. -/
theorem depolarizingChoi_eq {p : ℝ} (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    depolarizingChoi p = rankOneSum ![(1 - p) / 2, p / 4, p / 4, p / 4, p / 4]
      ![omega, ket 0 0, ket 0 1, ket 1 0, ket 1 1] := by
  have hs : (Real.sqrt (1 - 3 * p / 4) : ℂ) ^ 2 = 1 - 3 * p / 4 := by
    rw [← Complex.ofReal_pow, Real.sq_sqrt (by linarith)]; push_cast; ring
  have ht : (Real.sqrt p : ℂ) ^ 2 = p := by
    rw [← Complex.ofReal_pow, Real.sq_sqrt hp0]
  have h4 : Real.sqrt 4 = 2 := by
    rw [show (4 : ℝ) = 2 ^ 2 by norm_num, Real.sqrt_sq (by norm_num)]
  ext ⟨a, b⟩ ⟨a', b'⟩
  fin_cases a <;> fin_cases b <;> fin_cases a' <;> fin_cases b' <;>
    simp [depolarizingChoi, krausTerm, depolarizingKraus, phi, rankOneSum, vec, omega, ket,
      mul_apply, kroneckerMap_apply, conjTranspose_apply, vecMulVec_apply, Fin.sum_univ_succ,
      Fintype.sum_prod_type, one_apply, h4, map_ofNat] <;>
    (generalize (Real.sqrt (1 - 3 * p / 4) : ℂ) = S at hs ⊢
     generalize (Real.sqrt p : ℂ) = T at ht ⊢
     ring_nf; try simp only [hs, ht, Complex.I_sq]; try ring_nf)

/-- `I_max(D_p) = log2[(1 - p) d² + p]`, Eq. (S199): the optimum of Eq. (9) is `4 (1 - p) + p`.
Attained at `Y = (4 - 3p) I / 2`; certified by `Z = 4 Φ`. -/
theorem depolarizing_maxInformation {p : ℝ} (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    IsLeast (maxValues (depolarizingChoi p)) (4 * (1 - p) + p) := by
  rw [depolarizingChoi_eq hp0 hp1]
  refine isLeast_maxValues ((((4 - 3 * p) / 2 : ℝ) : ℂ) • 1)
    (PosSemidef.one.smul (Complex.zero_le_real.mpr (by linarith))) ?_ ?_ (c := ![2]) ![omega]
    (by intro i; fin_cases i; norm_num) ?_ ?_
  · have : maxMixed ⊗ₖ ((((4 - 3 * p) / 2 : ℝ) : ℂ) • (1 : Matrix (Fin 2) (Fin 2) ℂ)) -
        rankOneSum ![(1 - p) / 2, p / 4, p / 4, p / 4, p / 4]
          ![omega, ket 0 0, ket 0 1, ket 1 0, ket 1 1] =
        rankOneSum ![(1 - p) / 2, (1 - p) / 2, (1 - p) / 2] ![phiMinus, psiPlus, psiMinus] := by
      ext ⟨a, b⟩ ⟨a', b'⟩
      fin_cases a <;> fin_cases b <;> fin_cases a' <;> fin_cases b' <;> pctc_simp <;> ring
    rw [this]
    exact posSemidef_rankOneSum _ (by intro i; fin_cases i <;> simp <;> linarith)
  · pctc_simp; ring
  · ext b b'; fin_cases b <;> fin_cases b' <;> pctc_simp
  · pctc_simp; ring

/-- `I_doe(D_p) = -log2 p`, Eq. (S198): the optimum of Eq. (10) is `p`. Attained at `W = p I / 2`;
certified by `Z = 4 (I - Φ) / 3`. -/
theorem depolarizing_doeblinInformation {p : ℝ} (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    IsGreatest (doeValues (depolarizingChoi p)) p := by
  rw [depolarizingChoi_eq hp0 hp1]
  refine isGreatest_doeValues (((p / 2 : ℝ) : ℂ) • 1) ?_ ?_ ?_ (c := ![2 / 3, 2 / 3, 2 / 3])
    ![phiMinus, psiPlus, psiMinus] (by intro i; fin_cases i <;> norm_num) ?_ ?_
  · exact isHermitian_real_smul_one _
  · have : rankOneSum ![(1 - p) / 2, p / 4, p / 4, p / 4, p / 4]
          ![omega, ket 0 0, ket 0 1, ket 1 0, ket 1 1] -
        maxMixed ⊗ₖ (((p / 2 : ℝ) : ℂ) • (1 : Matrix (Fin 2) (Fin 2) ℂ)) =
        rankOneSum ![(1 - p) / 2] ![omega] := by
      ext ⟨a, b⟩ ⟨a', b'⟩
      fin_cases a <;> fin_cases b <;> fin_cases a' <;> fin_cases b' <;> pctc_simp <;> ring
    rw [this]
    exact posSemidef_rankOneSum _ (by intro i; fin_cases i; simp; linarith)
  · pctc_simp; ring
  · ext b b'; fin_cases b <;> fin_cases b' <;> pctc_simp <;> ring
  · pctc_simp; ring

/-! ## The erasure channel -/

/-- `J = (1 - p) Φ_q + p π ⊗ |e⟩⟨e|`. -/
theorem erasureChoi_eq {p : ℝ} (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    erasureChoi p = rankOneSum ![(1 - p) / 2, p / 2, p / 2] ![omegaQ, flag0, flag1] := by
  have hs : (Real.sqrt (1 - p) : ℂ) ^ 2 = 1 - p := by
    rw [← Complex.ofReal_pow, Real.sq_sqrt (by linarith)]; push_cast; ring
  have ht : (Real.sqrt p : ℂ) ^ 2 = p := by
    rw [← Complex.ofReal_pow, Real.sq_sqrt hp0]
  ext ⟨a, b⟩ ⟨a', b'⟩
  fin_cases a <;> fin_cases b <;> fin_cases a' <;> fin_cases b' <;>
    simp [erasureChoi, krausTerm, erasureKraus, phi, rankOneSum, vec, omegaQ, flag0, flag1,
      mul_apply, kroneckerMap_apply, conjTranspose_apply, vecMulVec_apply, Fin.sum_univ_succ,
      Fintype.sum_prod_type, one_apply] <;>
    (generalize (Real.sqrt (1 - p) : ℂ) = S at hs ⊢
     generalize (Real.sqrt p : ℂ) = T at ht ⊢
     ring_nf; try simp only [hs, ht]; try ring_nf)

/-- `I_max(E_p) = log2[(1 - p) d² + p]`, Eq. (S216): the optimum of Eq. (9) is `4 (1 - p) + p`.
Attained at `Y = 2 (1 - p) P_q + p |e⟩⟨e|`; certified by `Z = 4 Φ_q + I ⊗ |e⟩⟨e|`. -/
theorem erasure_maxInformation {p : ℝ} (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    IsLeast (maxValues (erasureChoi p)) (4 * (1 - p) + p) := by
  rw [erasureChoi_eq hp0 hp1]
  refine isLeast_maxValues (realDiagonal ![2 * (1 - p), 2 * (1 - p), p])
    (posSemidef_realDiagonal (by intro i; fin_cases i <;> simp <;> linarith)) ?_ ?_
    (c := ![2, 1, 1]) ![omegaQ, flag0, flag1] (by intro i; fin_cases i <;> norm_num) ?_ ?_
  · have : maxMixed ⊗ₖ realDiagonal ![2 * (1 - p), 2 * (1 - p), p] -
        rankOneSum ![(1 - p) / 2, p / 2, p / 2] ![omegaQ, flag0, flag1] =
        rankOneSum ![(1 - p) / 2, (1 - p) / 2, (1 - p) / 2] ![phiMinusQ, psiPlusQ, psiMinusQ] := by
      ext ⟨a, b⟩ ⟨a', b'⟩
      fin_cases a <;> fin_cases b <;> fin_cases a' <;> fin_cases b' <;> pctc_simp <;> pctc_finish
    rw [this]
    exact posSemidef_rankOneSum _ (by intro i; fin_cases i <;> simp <;> linarith)
  · pctc_simp <;> pctc_finish
  · ext b b'; fin_cases b <;> fin_cases b' <;> pctc_simp <;> pctc_finish
  · pctc_simp <;> pctc_finish

/-- `I_doe(E_p) = -log2 p`, Eq. (S215): the optimum of Eq. (10) is `p`. Attained at
`W = p |e⟩⟨e|`; certified by `Z = 4 (I ⊗ P_q - Φ_q) / 3 + I ⊗ |e⟩⟨e|`. -/
theorem erasure_doeblinInformation {p : ℝ} (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    IsGreatest (doeValues (erasureChoi p)) p := by
  rw [erasureChoi_eq hp0 hp1]
  refine isGreatest_doeValues (realDiagonal ![0, 0, p]) (isHermitian_realDiagonal _) ?_ ?_
    (c := ![2 / 3, 2 / 3, 2 / 3, 1, 1]) ![phiMinusQ, psiPlusQ, psiMinusQ, flag0, flag1]
    (by intro i; fin_cases i <;> norm_num) ?_ ?_
  · have : rankOneSum ![(1 - p) / 2, p / 2, p / 2] ![omegaQ, flag0, flag1] -
        maxMixed ⊗ₖ realDiagonal ![0, 0, p] = rankOneSum ![(1 - p) / 2] ![omegaQ] := by
      ext ⟨a, b⟩ ⟨a', b'⟩
      fin_cases a <;> fin_cases b <;> fin_cases a' <;> fin_cases b' <;> pctc_simp <;> pctc_finish
    rw [this]
    exact posSemidef_rankOneSum _ (by intro i; fin_cases i; simp; linarith)
  · pctc_simp <;> pctc_finish
  · ext b b'; fin_cases b <;> fin_cases b' <;> pctc_simp <;> pctc_finish
  · pctc_simp <;> pctc_finish

/-! ## The amplitude-damping channel -/

/-- The Kraus operators of `Channel.amplitudeDamping(γ)`, in Scotty's order. -/
def amplitudeDampingKraus (γ : ℝ) : Fin 2 → Matrix (Fin 2) (Fin 2) ℂ
  | 0 => !![1, 0; 0, (Real.sqrt (1 - γ) : ℂ)]
  | 1 => !![0, (Real.sqrt γ : ℂ); 0, 0]

/-- `Capacities.choiState(Channel.amplitudeDamping(γ))`. -/
def amplitudeDampingChoi (γ : ℝ) : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ :=
  ∑ k, krausTerm (amplitudeDampingKraus γ k)

/-- `J = (|v⟩⟨v| + γ |10⟩⟨10|) / 2` with `|v⟩ = |00⟩ + √(1 - γ) |11⟩`. -/
theorem amplitudeDampingChoi_eq {γ : ℝ} (hγ0 : 0 ≤ γ) (hγ1 : γ ≤ 1) :
    amplitudeDampingChoi γ = rankOneSum ![1 / 2, γ / 2]
      ![vec !![1, 0; 0, (Real.sqrt (1 - γ) : ℂ)], ket 1 0] := by
  have hs : (Real.sqrt γ : ℂ) ^ 2 = γ := by
    rw [← Complex.ofReal_pow, Real.sq_sqrt hγ0]
  ext ⟨a, b⟩ ⟨a', b'⟩
  fin_cases a <;> fin_cases b <;> fin_cases a' <;> fin_cases b' <;>
    simp [amplitudeDampingChoi, krausTerm, amplitudeDampingKraus, phi, rankOneSum, vec, ket,
      mul_apply, kroneckerMap_apply, conjTranspose_apply, vecMulVec_apply, Fin.sum_univ_succ,
      Fintype.sum_prod_type, one_apply] <;>
    (generalize (Real.sqrt γ : ℂ) = T at hs ⊢
     ring_nf; try simp only [hs]; try ring_nf)

/-- `I_max(A_γ) = log2[2 (1 + √(1 - γ)) - γ]`, Eq. (S223): the optimum of Eq. (9) is
`2 (1 + k) - γ = (1 + k)²`, `k = √(1 - γ)`. Attained at `Y = diag(1 + k, k² + k)`; certified by
`Z = 4 Φ`. -/
theorem amplitudeDamping_maxInformation {γ : ℝ} (hγ0 : 0 ≤ γ) (hγ1 : γ ≤ 1) :
    IsLeast (maxValues (amplitudeDampingChoi γ)) (2 * (1 + Real.sqrt (1 - γ)) - γ) := by
  rw [amplitudeDampingChoi_eq hγ0 hγ1]
  have hk0 : 0 ≤ Real.sqrt (1 - γ) := Real.sqrt_nonneg _
  have hk : (Real.sqrt (1 - γ) : ℂ) ^ 2 = 1 - γ := by
    rw [← Complex.ofReal_pow, Real.sq_sqrt (by linarith)]; push_cast; ring
  refine isLeast_maxValues
    (realDiagonal ![1 + Real.sqrt (1 - γ), Real.sqrt (1 - γ) ^ 2 + Real.sqrt (1 - γ)])
    (posSemidef_realDiagonal (by intro i; fin_cases i <;> simp <;> positivity)) ?_ ?_
    (c := ![2]) ![omega] (by intro i; fin_cases i; norm_num) ?_ ?_
  · have : maxMixed ⊗ₖ realDiagonal ![1 + Real.sqrt (1 - γ), Real.sqrt (1 - γ) ^ 2 + Real.sqrt (1 - γ)] -
        rankOneSum ![1 / 2, γ / 2] ![vec !![1, 0; 0, (Real.sqrt (1 - γ) : ℂ)], ket 1 0] =
        rankOneSum ![Real.sqrt (1 - γ) / 2, (Real.sqrt (1 - γ) + Real.sqrt (1 - γ) ^ 2) / 2,
          (Real.sqrt (1 - γ) ^ 2 + Real.sqrt (1 - γ)) / 2] ![phiMinus, ket 1 0, ket 0 1] := by
      ext ⟨a, b⟩ ⟨a', b'⟩
      fin_cases a <;> fin_cases b <;> fin_cases a' <;> fin_cases b' <;> pctc_simp <;>
        (generalize (Real.sqrt (1 - γ) : ℂ) = K at hk ⊢
         ring_nf; try simp only [hk]; try ring_nf)
    rw [this]
    exact posSemidef_rankOneSum _ (by intro i; fin_cases i <;> simp <;> positivity)
  · pctc_simp
    generalize (Real.sqrt (1 - γ) : ℂ) = K at hk ⊢
    linear_combination hk
  · ext b b'; fin_cases b <;> fin_cases b' <;> pctc_simp <;> pctc_finish
  · pctc_simp
    generalize (Real.sqrt (1 - γ) : ℂ) = K at hk ⊢
    linear_combination hk

/-- `I_doe(A_γ) = -2 log2[1 - √(1 - γ)]`, Eq. (S222): the optimum of Eq. (10) is `(1 - k)²`,
`k = √(1 - γ)`. Attained at `W = diag(1 - k, k² - k)`; certified by `Z = 2 |u⟩⟨u|`,
`|u⟩ = |00⟩ - |11⟩`. -/
theorem amplitudeDamping_doeblinInformation {γ : ℝ} (hγ0 : 0 ≤ γ) (hγ1 : γ ≤ 1) :
    IsGreatest (doeValues (amplitudeDampingChoi γ)) ((1 - Real.sqrt (1 - γ)) ^ 2) := by
  rw [amplitudeDampingChoi_eq hγ0 hγ1]
  have hk0 : 0 ≤ Real.sqrt (1 - γ) := Real.sqrt_nonneg _
  have hk1 : Real.sqrt (1 - γ) ≤ 1 := Real.sqrt_le_one.mpr (by linarith)
  have hk : (Real.sqrt (1 - γ) : ℂ) ^ 2 = 1 - γ := by
    rw [← Complex.ofReal_pow, Real.sq_sqrt (by linarith)]; push_cast; ring
  refine isGreatest_doeValues
    (realDiagonal ![1 - Real.sqrt (1 - γ), Real.sqrt (1 - γ) ^ 2 - Real.sqrt (1 - γ)])
    (isHermitian_realDiagonal _) ?_ ?_ (c := ![2]) ![phiMinus] (by intro i; fin_cases i; norm_num)
    ?_ ?_
  · have : rankOneSum ![1 / 2, γ / 2] ![vec !![1, 0; 0, (Real.sqrt (1 - γ) : ℂ)], ket 1 0] -
        maxMixed ⊗ₖ realDiagonal ![1 - Real.sqrt (1 - γ), Real.sqrt (1 - γ) ^ 2 - Real.sqrt (1 - γ)] =
        rankOneSum ![Real.sqrt (1 - γ) / 2, (Real.sqrt (1 - γ) - Real.sqrt (1 - γ) ^ 2) / 2,
          (Real.sqrt (1 - γ) - Real.sqrt (1 - γ) ^ 2) / 2] ![omega, ket 1 0, ket 0 1] := by
      ext ⟨a, b⟩ ⟨a', b'⟩
      fin_cases a <;> fin_cases b <;> fin_cases a' <;> fin_cases b' <;> pctc_simp <;>
        (generalize (Real.sqrt (1 - γ) : ℂ) = K at hk ⊢
         ring_nf; try simp only [hk]; try ring_nf)
    rw [this]
    refine posSemidef_rankOneSum _ (by intro i; fin_cases i <;> simp <;> nlinarith)
  · pctc_simp <;> pctc_finish
  · ext b b'; fin_cases b <;> fin_cases b' <;> pctc_simp <;> pctc_finish
  · pctc_simp <;> pctc_finish

end

end Pctc
