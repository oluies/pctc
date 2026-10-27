import Pctc.Basic

/-!
# The noiseless P-CTC loop reduces to the Lloyd et al. operator

`PctcLoop.evaluate` in `src/main/scala/pctc/PctcLoop.scala` computes Eq. (1) of arXiv:2509.08965v3.
Its scaladoc for `loopOperator` claims that with a noiseless loop and a unitary bipartite map the
result is `rho ↦ C rho Cᴴ / Tr[C rho Cᴴ]` with `C = Tr_CTC[U]` (arXiv:1005.2219), at a
post-selection probability of `Tr[C rho Cᴴ] / d²`. The test suite checks that numerically; this file
proves it.

`loopNumerator` follows `evaluate` step by step, with the reference system `R` trivial,
`N = id` and `T(X) = U X Uᴴ`, on `E = F = S` (chronology respecting) and `B = A = K` (the CTC).
`U` does not have to be unitary: the identity is purely algebraic.

The headline results are `loopNumerator_eq` and `loop_normalized_eq`.
-/

namespace Pctc

open Matrix
open scoped Kronecker

variable {S K : Type*} [Fintype S] [DecidableEq S] [Fintype K] [DecidableEq K]

noncomputable section

/-- The Lloyd et al. loop operator `C = Tr_CTC[U]`, `PctcLoop.loopOperator`. -/
def loopOperator (U : Matrix (S × K) (S × K) ℂ) : Matrix S S ℂ := traceRight U

/-! ## The steps of `PctcLoop.evaluate` -/

/-- `rhoRE.tensor(loopResource)`, on `(E, B, A')`. With `N = id` the loop resource is `Φ`. -/
def loopInput (ρ : Matrix S S ℂ) : Matrix ((S × K) × K) ((S × K) × K) ℂ :=
  reindex (Equiv.prodAssoc S K K).symm (Equiv.prodAssoc S K K).symm (ρ ⊗ₖ phi)

/-- `t.applyTo(input, Seq(1, 2))`, on `(F, A, A')`. -/
def loopEvolved (U : Matrix (S × K) (S × K) ℂ) (ρ : Matrix S S ℂ) :
    Matrix ((S × K) × K) ((S × K) × K) ℂ :=
  (U ⊗ₖ (1 : Matrix K K ℂ)) * loopInput ρ * (U ⊗ₖ (1 : Matrix K K ℂ))ᴴ

/-- `projector.multiply(evolved.matrix)`, with `projector = I ⊗ Φ_{AA'}`. -/
def loopWeighted (U : Matrix (S × K) (S × K) ℂ) (ρ : Matrix S S ℂ) :
    Matrix (S × (K × K)) (S × (K × K)) ℂ :=
  ((1 : Matrix S S ℂ) ⊗ₖ phi) *
    reindex (Equiv.prodAssoc S K K) (Equiv.prodAssoc S K K) (loopEvolved U ρ)

/-- `weighted.partialTrace(Seq(2, 3))`: the numerator of Eq. (1). -/
def loopNumerator (U : Matrix (S × K) (S × K) ℂ) (ρ : Matrix S S ℂ) : Matrix S S ℂ :=
  traceRight (loopWeighted U ρ)

/-! ## Entries -/

omit [Fintype S] [DecidableEq S] in
theorem loopInput_apply (ρ : Matrix S S ℂ) (s t : S) (b c a e : K) :
    loopInput (K := K) ρ ((s, b), a) ((t, c), e) = ρ s t * phi (b, a) (c, e) := rfl

omit [DecidableEq S] in
theorem loopEvolved_apply (U : Matrix (S × K) (S × K) ℂ) (ρ : Matrix S S ℂ)
    (i j : S) (b b' a a' : K) :
    loopEvolved U ρ ((i, b), b') ((j, a), a') =
      (Fintype.card K : ℂ)⁻¹ * ∑ s, ∑ t, U (i, b) (s, b') * ρ s t * star (U (j, a) (t, a')) := by
  simp [loopEvolved, mul_apply, Fintype.sum_prod_type, kroneckerMap_apply, one_apply,
    loopInput_apply, phi, conjTranspose_apply, apply_ite (starRingEnd ℂ)]
  simp_rw [Finset.sum_mul, Finset.mul_sum]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun _ _ => Finset.sum_congr rfl fun _ _ => ?_
  ring

/-! ## The reduction -/

/-- The noiseless loop is `C ρ Cᴴ`, scaled by `1 / d²`. -/
theorem loopNumerator_eq (U : Matrix (S × K) (S × K) ℂ) (ρ : Matrix S S ℂ) :
    loopNumerator U ρ =
      ((Fintype.card K : ℂ)⁻¹ ^ 2) • (loopOperator U * ρ * (loopOperator U)ᴴ) := by
  ext i j
  simp [loopNumerator, traceRight, loopWeighted, mul_apply, Fintype.sum_prod_type,
    kroneckerMap_apply, one_apply, phi, loopEvolved_apply, loopOperator, conjTranspose_apply]
  simp_rw [← Finset.mul_sum]
  rw [← mul_assoc, ← mul_inv, ← sq]
  congr 1
  simp_rw [Finset.sum_mul, Finset.mul_sum]
  -- both sides are one sum over (a, b, s, t), in the orders (a, b, s, t) and (t, s, b, a)
  conv_lhs => rw [Finset.sum_comm]
  conv_lhs => enter [2, b]; rw [Finset.sum_comm]
  conv_lhs => enter [2, b, 2, s]; rw [Finset.sum_comm]
  conv_lhs => rw [Finset.sum_comm]
  conv_lhs => enter [2, s]; rw [Finset.sum_comm]
  conv_lhs => rw [Finset.sum_comm]

/-- The denominator of Eq. (1), `numerator.realTrace`: `Tr[C ρ Cᴴ] / d²`. -/
theorem trace_loopNumerator (U : Matrix (S × K) (S × K) ℂ) (ρ : Matrix S S ℂ) :
    (loopNumerator U ρ).trace =
      (Fintype.card K : ℂ)⁻¹ ^ 2 * (loopOperator U * ρ * (loopOperator U)ᴴ).trace := by
  rw [loopNumerator_eq, trace_smul, smul_eq_mul]

/-- After the loop divides out its post-selection probability, the output is the Lloyd et al.
evolution `ρ ↦ C ρ Cᴴ / Tr[C ρ Cᴴ]`. -/
theorem loop_normalized_eq [Nonempty K] (U : Matrix (S × K) (S × K) ℂ) (ρ : Matrix S S ℂ) :
    (loopNumerator U ρ).trace⁻¹ • loopNumerator U ρ =
      (loopOperator U * ρ * (loopOperator U)ᴴ).trace⁻¹ •
        (loopOperator U * ρ * (loopOperator U)ᴴ) := by
  have hd : (Fintype.card K : ℂ)⁻¹ ^ 2 ≠ 0 := by simp
  rw [trace_loopNumerator, loopNumerator_eq, smul_smul, mul_inv, mul_right_comm,
    inv_mul_cancel₀ hd, one_mul]

/-! ## Demonstration (a): the grandfather paradox -/

/-- `RY(θ)`. -/
def ry (θ : ℝ) : Matrix (Fin 2) (Fin 2) ℂ :=
  !![(Real.cos (θ / 2) : ℂ), -(Real.sin (θ / 2) : ℂ); (Real.sin (θ / 2) : ℂ), (Real.cos (θ / 2) : ℂ)]

/-- `RY(θ)` on the CTC qubit of an otherwise empty loop: no chronology-respecting system. -/
def grandfatherUnitary (θ : ℝ) : Matrix (Unit × Fin 2) (Unit × Fin 2) ℂ :=
  fun p q => ry θ p.2 q.2

/-- The post-selection probability in demonstration (a) is `cos²(θ/2)`. -/
theorem grandfather_postSelection (θ : ℝ) :
    (loopNumerator (grandfatherUnitary θ) 1).trace = ((Real.cos (θ / 2) ^ 2 : ℝ) : ℂ) := by
  rw [trace_loopNumerator]
  simp [loopOperator, traceRight, grandfatherUnitary, ry, Matrix.trace, Fin.sum_univ_two,
    mul_apply, conjTranspose_apply, -Complex.ofReal_cos, -Complex.ofReal_sin]
  ring

/-- At `θ = π` the rotation is a bit flip and the branch has probability exactly zero. -/
theorem grandfather_postSelection_pi :
    (loopNumerator (grandfatherUnitary Real.pi) 1).trace = 0 := by
  simp [grandfather_postSelection]

end

end Pctc
