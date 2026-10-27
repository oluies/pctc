import Mathlib

/-!
# Definitions shared by the pctc proofs

The partial trace and the maximally entangled state, as Scotty's `DensityMatrix` builds them.
-/

namespace Pctc

open Matrix

variable {S K : Type*} [Fintype S] [DecidableEq S] [Fintype K] [DecidableEq K]

noncomputable section

/-- The partial trace over the right tensor factor, `Tr_K`. -/
def traceRight (M : Matrix (S × K) (S × K) ℂ) : Matrix S S ℂ :=
  fun i j => ∑ k, M (i, k) (j, k)

/-- The standard maximally entangled state `Φ = |Ω⟩⟨Ω| / d` on `K × K`, `|Ω⟩ = ∑ₖ |k k⟩`.
`DensityMatrix.maximallyEntangled` in Scotty. -/
def phi : Matrix (K × K) (K × K) ℂ :=
  fun p q => (Fintype.card K : ℂ)⁻¹ * (if p.1 = p.2 then 1 else 0) * (if q.1 = q.2 then 1 else 0)

end

end Pctc
