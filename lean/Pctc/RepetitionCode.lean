import Mathlib

/-!
# The repetition code, decoded by majority vote

A side section, not part of arXiv:2509.08965. Cat qubits barely ever bit-flip, so Alice & Bob protect
them against phase flips alone, with a repetition code (IACR ePrint 2026/2665). This file proves the
two facts the page's simulation rests on:

* majority vote on `2t + 1` copies corrects every pattern of at most `t` flips, and no more;
* with each copy flipped independently with probability `p`, the decoder fails with exactly the
  binomial tail probability, `RepetitionCode.phaseFlipFailure` in the Scala.

An error pattern is the set of positions that flipped.
-/

namespace Pctc

open Finset

variable {n : ℕ}

/-- The word read back when codeword `b, …, b` suffers the flips in `e`. -/
def received (b : Bool) (e : Finset (Fin n)) : Fin n → Bool :=
  fun i => if i ∈ e then !b else b

/-- Majority vote: `true` when more than half of the `n` bits are. -/
def majority (r : Fin n → Bool) : Bool :=
  decide (n < 2 * #(univ.filter fun i => r i = true))

theorem majority_received_false (e : Finset (Fin n)) :
    majority (received false e) = decide (n < 2 * #e) := by
  have : (univ.filter fun i => received false e i = true) = e := by
    ext i; simp [received]
  simp [majority, this]

theorem majority_received_true (e : Finset (Fin n)) :
    majority (received true e) = decide (n < 2 * (n - #e)) := by
  have : (univ.filter fun i => received true e i = true) = eᶜ := by
    ext i; simp [received]
  simp [majority, this, card_compl]

/-- Majority vote on `2t + 1` copies recovers the codeword from any `t` or fewer flips. -/
theorem majority_corrects {t : ℕ} (b : Bool) (e : Finset (Fin (2 * t + 1))) (he : #e ≤ t) :
    majority (received b e) = b := by
  cases b
  · rw [majority_received_false]; simp; omega
  · rw [majority_received_true]; simp; omega

/-- … and `t + 1` flips in the wrong places defeat it, so `t` is the most it corrects. -/
theorem majority_fails {t : ℕ} (b : Bool) :
    ∃ e : Finset (Fin (2 * t + 1)), #e = t + 1 ∧ majority (received b e) ≠ b := by
  obtain ⟨e, -, he⟩ := exists_subset_card_eq (s := (univ : Finset (Fin (2 * t + 1))))
    (n := t + 1) (by simp; omega)
  refine ⟨e, he, ?_⟩
  cases b
  · rw [majority_received_false]; simp; omega
  · rw [majority_received_true]; simp; omega

/-- The probability that majority vote misreads codeword `false`, when each of the `n` copies
flips independently with probability `p`: the total weight of the patterns it misreads. -/
noncomputable def failureProbability (n : ℕ) (p : ℝ) : ℝ :=
  ∑ e ∈ (univ : Finset (Fin n)).powerset,
    if majority (received false e) then p ^ #e * (1 - p) ^ (n - #e) else 0

/-- Grouping the patterns by weight gives the binomial tail: more than half the copies flipped. -/
theorem failureProbability_eq (n : ℕ) (p : ℝ) :
    failureProbability n p =
      ∑ m ∈ range (n + 1), if n < 2 * m then (n.choose m : ℝ) * p ^ m * (1 - p) ^ (n - m) else 0 := by
  unfold failureProbability
  simp_rw [majority_received_false]
  rw [sum_powerset_apply_card (fun m => if decide (n < 2 * m) then p ^ m * (1 - p) ^ (n - m) else 0)]
  simp only [card_univ, Fintype.card_fin, nsmul_eq_mul]
  refine sum_congr rfl fun m _ => ?_
  split_ifs <;> simp_all [mul_assoc]

end Pctc
