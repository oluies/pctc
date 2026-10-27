import Pctc.Basic

/-!
# For the erasure channel the optimum can be taken block diagonal

`Bloch.stateFamily` and `Bloch.hermitianFamily` sweep only block-diagonal qutrits for the erasure
channel: a qubit block on `|0⟩, |1⟩` plus a weight on the flag `|e⟩ = |2⟩`, with no coherence
between the two. The scaladoc of `Bloch.blockQutrit` argues that this loses nothing. This file
proves it.

With `λ` folded into the variable, the two programs of `Capacities` are

* `I_max`: minimise `Tr Y` over `Y ≥ 0` with `π ⊗ Y - J ≥ 0`, Eq. (9);
* `I_doe`: maximise `Tr W` over Hermitian `W` with `J - π ⊗ W ≥ 0`, Eq. (10);

where `J` is the Choi operator on `(A', B)`. Let `D` be `+1` on the qubit block and `-1` on the
flag, and `pinch Y = (Y + D Y D) / 2`, which keeps the two diagonal blocks of `Y` and zeroes the
coherences between them. If `(I ⊗ D) J (I ⊗ D) = J`, then pinching preserves both constraints and
the trace (`maxInformation_pinch`, `doeblinInformation_pinch`), so every value either program
reaches it also reaches at a block-diagonal operator.

The erasure channel's Kraus operators each map into a single block, `D K = ±K`, which is what makes
its Choi operator invariant (`erasureChoi_invariant`). The scaladoc averages over the whole phase
group `I₂ + e^{iφ}|e⟩⟨e|`; the two-element subgroup `φ ∈ {0, π}` is already enough, and is what is
used here.

The headline results are `erasure_maxInformation_blockDiagonal` and
`erasure_doeblinInformation_blockDiagonal`. Neither depends on `π`, which is `I / 2` in
`Capacities`, nor on the erasure probability `p`.
-/

namespace Pctc

open Matrix
open scoped Kronecker ComplexOrder

noncomputable section

variable {n : Type*} [Fintype n]

/-! ## Averaging with an involution -/

/-- `(M + G M G) / 2`. With `G` a Hermitian involution this averages `M` over the group `{1, G}`. -/
def average (G M : Matrix n n ℂ) : Matrix n n ℂ := ((1 / 2 : ℝ) : ℂ) • (M + G * M * G)

theorem average_sub (G M N : Matrix n n ℂ) : average G (M - N) = average G M - average G N := by
  simp only [average, mul_sub, sub_mul]
  rw [← smul_sub]
  congr 1
  abel

theorem average_eq_self {G N : Matrix n n ℂ} (h : G * N * G = N) : average G N = N := by
  rw [average, h, ← two_smul ℂ N, smul_smul]
  norm_num

theorem posSemidef_average {G M : Matrix n n ℂ} (hG : Gᴴ = G) (hM : M.PosSemidef) :
    (average G M).PosSemidef := by
  have hGMG : (G * M * G).PosSemidef := by
    simpa [hG] using hM.mul_mul_conjTranspose_same G
  exact (hM.add hGMG).smul (Complex.zero_le_real.mpr (by norm_num))

theorem isHermitian_average {G M : Matrix n n ℂ} (hG : Gᴴ = G) (hM : M.IsHermitian) :
    (average G M).IsHermitian := by
  unfold Matrix.IsHermitian average
  rw [conjTranspose_smul, conjTranspose_add, conjTranspose_mul, conjTranspose_mul, hG, hM.eq,
    Complex.star_def, Complex.conj_ofReal, mul_assoc]

/-! ## The block sign and the pinching it induces -/

variable {B : Type*} [Fintype B] [DecidableEq B] (flag : B → Prop) [DecidablePred flag]

/-- `D`: `+1` off the flag block, `-1` on it. -/
def blockSign : Matrix B B ℂ := diagonal fun b => if flag b then -1 else 1

omit [Fintype B] in
theorem blockSign_conjTranspose : (blockSign flag)ᴴ = blockSign flag := by
  rw [blockSign, diagonal_conjTranspose]
  congr 1
  ext b
  by_cases h : flag b <;> simp [h]

/-- `Y` has no entries between the flag block and its complement. -/
def IsBlockDiag (Y : Matrix B B ℂ) : Prop := ∀ i j, ¬(flag i ↔ flag j) → Y i j = 0

/-- `(Y + D Y D) / 2`, the block-diagonal part of `Y`. -/
def pinch (Y : Matrix B B ℂ) : Matrix B B ℂ := average (blockSign flag) Y

theorem pinch_apply (Y : Matrix B B ℂ) (i j : B) :
    pinch flag Y i j = if (flag i ↔ flag j) then Y i j else 0 := by
  by_cases hi : flag i <;> by_cases hj : flag j <;>
    simp [pinch, average, blockSign, mul_diagonal, diagonal_mul, hi, hj] <;> ring

theorem pinch_isBlockDiag (Y : Matrix B B ℂ) : IsBlockDiag flag (pinch flag Y) := by
  intro i j h
  rw [pinch_apply, ite_eq_right_iff]
  exact fun h' => absurd h' h

theorem trace_pinch (Y : Matrix B B ℂ) : (pinch flag Y).trace = Y.trace := by
  simp only [Matrix.trace, diag_apply, pinch_apply, iff_self, ite_true]

/-! ## Pinching on the output factor -/

variable {A : Type*} [Fintype A] [DecidableEq A]

omit [Fintype A] [Fintype B] in
theorem kronecker_blockSign_conjTranspose :
    ((1 : Matrix A A ℂ) ⊗ₖ blockSign flag)ᴴ = (1 : Matrix A A ℂ) ⊗ₖ blockSign flag := by
  rw [conjTranspose_kronecker, conjTranspose_one, blockSign_conjTranspose]

theorem kronecker_pinch (π : Matrix A A ℂ) (Y : Matrix B B ℂ) :
    π ⊗ₖ pinch flag Y = average ((1 : Matrix A A ℂ) ⊗ₖ blockSign flag) (π ⊗ₖ Y) := by
  rw [pinch, average, average, kronecker_smul, kronecker_add, ← mul_kronecker_mul,
    ← mul_kronecker_mul, one_mul, mul_one]

/-- `I_max`: pinching keeps `Y ≥ 0`, `π ⊗ Y - J ≥ 0` and `Tr Y`. -/
theorem maxInformation_pinch (π : Matrix A A ℂ) (J : Matrix (A × B) (A × B) ℂ)
    (hJ : ((1 : Matrix A A ℂ) ⊗ₖ blockSign flag) * J * ((1 : Matrix A A ℂ) ⊗ₖ blockSign flag) = J)
    {Y : Matrix B B ℂ} (hY : Y.PosSemidef) (h : (π ⊗ₖ Y - J).PosSemidef) :
    (pinch flag Y).PosSemidef ∧ (π ⊗ₖ pinch flag Y - J).PosSemidef ∧
      (pinch flag Y).trace = Y.trace := by
  refine ⟨posSemidef_average (blockSign_conjTranspose flag) hY, ?_, trace_pinch flag Y⟩
  rw [kronecker_pinch, ← average_eq_self hJ, ← average_sub]
  exact posSemidef_average (kronecker_blockSign_conjTranspose flag) h

/-- `I_doe`: pinching keeps `W` Hermitian, `J - π ⊗ W ≥ 0` and `Tr W`. -/
theorem doeblinInformation_pinch (π : Matrix A A ℂ) (J : Matrix (A × B) (A × B) ℂ)
    (hJ : ((1 : Matrix A A ℂ) ⊗ₖ blockSign flag) * J * ((1 : Matrix A A ℂ) ⊗ₖ blockSign flag) = J)
    {W : Matrix B B ℂ} (hW : W.IsHermitian) (h : (J - π ⊗ₖ W).PosSemidef) :
    (pinch flag W).IsHermitian ∧ (J - π ⊗ₖ pinch flag W).PosSemidef ∧
      (pinch flag W).trace = W.trace := by
  refine ⟨isHermitian_average (blockSign_conjTranspose flag) hW, ?_, trace_pinch flag W⟩
  rw [kronecker_pinch, ← average_eq_self hJ, ← average_sub]
  exact posSemidef_average (kronecker_blockSign_conjTranspose flag) h

/-! ## Choi operators of block-preserving channels -/

/-- If `G X = c X` with `|c| = 1`, conjugating `X P Xᴴ` by `G` leaves it unchanged. -/
theorem sandwich_eq_self {m p : Type*} [Fintype m] [Fintype p] {G : Matrix m m ℂ} {X : Matrix m p ℂ}
    {P : Matrix p p ℂ} {c : ℂ} (hl : G * X = c • X) (hr : Xᴴ * G = star c • Xᴴ)
    (hc : star c * c = 1) : G * (X * P * Xᴴ) * G = X * P * Xᴴ := by
  calc G * (X * P * Xᴴ) * G = G * X * P * (Xᴴ * G) := by simp only [Matrix.mul_assoc]
    _ = X * P * Xᴴ := by
      rw [hl, hr, Matrix.smul_mul, Matrix.smul_mul, Matrix.mul_smul, smul_smul, mul_comm, hc,
        one_smul]

/-- One Kraus branch of the Choi operator, `(I ⊗ K) Φ (I ⊗ K)ᴴ`. -/
def krausTerm (K : Matrix B A ℂ) : Matrix (A × B) (A × B) ℂ :=
  (((1 : Matrix A A ℂ) ⊗ₖ K) * (phi : Matrix (A × A) (A × A) ℂ) : Matrix (A × B) (A × A) ℂ) *
    conjTranspose ((1 : Matrix A A ℂ) ⊗ₖ K)

/-- A channel whose every Kraus operator maps into one block, `D K = ±K`, has a Choi operator
invariant under `I ⊗ D`. -/
theorem choi_invariant {ι : Type*} [Fintype ι] (K : ι → Matrix B A ℂ) (c : ι → ℂ)
    (hc : ∀ k, star (c k) * c k = 1) (hK : ∀ k, blockSign flag * K k = c k • K k) :
    ((1 : Matrix A A ℂ) ⊗ₖ blockSign flag) *
        (∑ k, krausTerm (K k)) * ((1 : Matrix A A ℂ) ⊗ₖ blockSign flag) =
      ∑ k, krausTerm (K k) := by
  rw [Finset.mul_sum, Finset.sum_mul]
  refine Finset.sum_congr rfl fun k _ => ?_
  unfold krausTerm
  have hl : ((1 : Matrix A A ℂ) ⊗ₖ blockSign flag) * ((1 : Matrix A A ℂ) ⊗ₖ K k) =
      c k • ((1 : Matrix A A ℂ) ⊗ₖ K k) := by
    rw [← mul_kronecker_mul, one_mul, hK, kronecker_smul]
  have hr : conjTranspose ((1 : Matrix A A ℂ) ⊗ₖ K k) * ((1 : Matrix A A ℂ) ⊗ₖ blockSign flag) =
      star (c k) • conjTranspose ((1 : Matrix A A ℂ) ⊗ₖ K k) := by
    rw [← conjTranspose_smul, ← hl, conjTranspose_mul, kronecker_blockSign_conjTranspose]
  exact sandwich_eq_self hl hr (hc k)

/-! ## The erasure channel -/

/-- The flag `|e⟩`: basis element 2 of the qutrit output, as in Scotty's `Channel.erasure`. -/
def erasureFlag (b : Fin 3) : Prop := b = 2

instance : DecidablePred erasureFlag := fun b => inferInstanceAs (Decidable (b = 2))

/-- The Kraus operators of `Channel.erasure(p)`, in Scotty's order. -/
def erasureKraus (p : ℝ) : Fin 3 → Matrix (Fin 3) (Fin 2) ℂ
  | 0 => (Real.sqrt (1 - p) : ℂ) • !![1, 0; 0, 1; 0, 0]
  | 1 => (Real.sqrt p : ℂ) • !![0, 0; 0, 0; 1, 0]
  | 2 => (Real.sqrt p : ℂ) • !![0, 0; 0, 0; 0, 1]

/-- `N(Φ_{A'A})` on `(A', B)`, `Capacities.choiState(Channel.erasure(p))`. -/
def erasureChoi (p : ℝ) : Matrix (Fin 2 × Fin 3) (Fin 2 × Fin 3) ℂ :=
  ∑ k, krausTerm (erasureKraus p k)

/-- The kept branch lies in the qubit block, the two erased branches on the flag. -/
def erasureSign : Fin 3 → ℂ
  | 0 => 1
  | 1 => -1
  | 2 => -1

theorem blockSign_mul_erasureKraus (p : ℝ) (k : Fin 3) :
    blockSign erasureFlag * erasureKraus p k = erasureSign k • erasureKraus p k := by
  ext i j
  fin_cases k <;> fin_cases i <;> fin_cases j <;>
    simp [blockSign, erasureFlag, erasureKraus, erasureSign, diagonal_mul]

/-- The kept branch: `(1 - p) Φ` on the qubit block. -/
theorem erasureChoi_kept {p : ℝ} (hp : p ≤ 1) (a b a' b' : Fin 2) :
    erasureChoi p (a, b.castSucc) (a', b'.castSucc) =
      ((1 - p) / 2 : ℝ) * (if a = b then 1 else 0) * (if a' = b' then 1 else 0) := by
  have h : (Real.sqrt (1 - p) : ℂ) * (Real.sqrt (1 - p) : ℂ) = ((1 - p : ℝ) : ℂ) := by
    rw [← Complex.ofReal_mul, Real.mul_self_sqrt (by linarith)]
  push_cast at h
  fin_cases a <;> fin_cases b <;> fin_cases a' <;> fin_cases b' <;>
    simp [erasureChoi, krausTerm, erasureKraus, phi, mul_apply, kroneckerMap_apply,
      conjTranspose_apply, Fin.sum_univ_three, Fin.sum_univ_two, Fintype.sum_prod_type, one_apply] <;>
    linear_combination (1 / 2 : ℂ) * h

/-- The erased branches: `p π ⊗ |e⟩⟨e|` on the flag. -/
theorem erasureChoi_flag {p : ℝ} (hp : 0 ≤ p) (a a' : Fin 2) :
    erasureChoi p (a, 2) (a', 2) = if a = a' then ((p / 2 : ℝ) : ℂ) else 0 := by
  have h : (Real.sqrt p : ℂ) * (Real.sqrt p : ℂ) = (p : ℂ) := by
    rw [← Complex.ofReal_mul, Real.mul_self_sqrt hp]
  fin_cases a <;> fin_cases a' <;>
    simp [erasureChoi, krausTerm, erasureKraus, phi, mul_apply, kroneckerMap_apply,
      conjTranspose_apply, Fin.sum_univ_three, Fin.sum_univ_two, Fintype.sum_prod_type, one_apply] <;>
    linear_combination (1 / 2 : ℂ) * h

theorem erasureChoi_invariant (p : ℝ) :
    ((1 : Matrix (Fin 2) (Fin 2) ℂ) ⊗ₖ blockSign erasureFlag) * erasureChoi p *
        ((1 : Matrix (Fin 2) (Fin 2) ℂ) ⊗ₖ blockSign erasureFlag) = erasureChoi p :=
  choi_invariant erasureFlag (erasureKraus p) erasureSign
    (fun k => by fin_cases k <;> simp [erasureSign]) (blockSign_mul_erasureKraus p)

/-- `I_max` of the erasure channel: every feasible `Y` has a block-diagonal feasible counterpart
with the same trace. -/
theorem erasure_maxInformation_blockDiagonal (p : ℝ) (π : Matrix (Fin 2) (Fin 2) ℂ)
    {Y : Matrix (Fin 3) (Fin 3) ℂ} (hY : Y.PosSemidef) (h : (π ⊗ₖ Y - erasureChoi p).PosSemidef) :
    ∃ Y', IsBlockDiag erasureFlag Y' ∧ Y'.PosSemidef ∧ (π ⊗ₖ Y' - erasureChoi p).PosSemidef ∧
      Y'.trace = Y.trace :=
  ⟨pinch erasureFlag Y, pinch_isBlockDiag erasureFlag Y,
    maxInformation_pinch erasureFlag π _ (erasureChoi_invariant p) hY h⟩

/-- `I_doe` of the erasure channel: every feasible `W` has a block-diagonal feasible counterpart
with the same trace. -/
theorem erasure_doeblinInformation_blockDiagonal (p : ℝ) (π : Matrix (Fin 2) (Fin 2) ℂ)
    {W : Matrix (Fin 3) (Fin 3) ℂ} (hW : W.IsHermitian) (h : (erasureChoi p - π ⊗ₖ W).PosSemidef) :
    ∃ W', IsBlockDiag erasureFlag W' ∧ W'.IsHermitian ∧ (erasureChoi p - π ⊗ₖ W').PosSemidef ∧
      W'.trace = W.trace :=
  ⟨pinch erasureFlag W, pinch_isBlockDiag erasureFlag W,
    doeblinInformation_pinch erasureFlag π _ (erasureChoi_invariant p) hW h⟩

end

end Pctc
