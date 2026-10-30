import Mathlib

/-!
# The repetition code, decoded by majority vote

A side section, not part of arXiv:2509.08965. Cat qubits barely ever bit-flip, so Alice & Bob protect
them against phase flips alone, with a repetition code (IACR ePrint 2026/2665). This file proves the
two facts the page's simulation rests on:

* majority vote on `2t + 1` copies corrects every pattern of at most `t` flips, and no more;
* with each copy flipped independently with probability `p`, the decoder fails with exactly the
  binomial tail probability, `RepetitionCode.phaseFlipFailureClosedForm` in the Scala.

It also proves what a cheaper decoder gives up. Recursive majority on nine copies, the majority of
three majorities of three, is what brokkr's `RecursiveMajority` builds: it corrects every pattern of
three flips but not every pattern of four, where exact majority on nine corrects four.

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

/-! ### Recursive majority on nine copies -/

/-- Majority of three bits. -/
def maj3 (a b c : Bool) : Bool := (a && b) || (b && c) || (a && c)

/-- The majority of the majorities of positions `0–2`, `3–5` and `6–8`. -/
def recursiveMajority (r : Fin 9 → Bool) : Bool :=
  maj3 (maj3 (r 0) (r 1) (r 2)) (maj3 (r 3) (r 4) (r 5)) (maj3 (r 6) (r 7) (r 8))

/-- Exact majority on nine copies corrects four flips: `majority_corrects` at `t = 4`. -/
theorem majority_nine_corrects (b : Bool) (e : Finset (Fin 9)) (he : #e ≤ 4) :
    majority (received b e) = b :=
  majority_corrects (t := 4) b e he

/-- The whole case analysis, over the codeword and the nine flip bits: 1024 cases, checked by the
kernel. -/
theorem recursiveMajority_bits : ∀ b e₀ e₁ e₂ e₃ e₄ e₅ e₆ e₇ e₈ : Bool,
    e₀.toNat + e₁.toNat + e₂.toNat + e₃.toNat + e₄.toNat + e₅.toNat + e₆.toNat + e₇.toNat + e₈.toNat ≤ 3 →
    maj3 (maj3 (b ^^ e₀) (b ^^ e₁) (b ^^ e₂)) (maj3 (b ^^ e₃) (b ^^ e₄) (b ^^ e₅))
      (maj3 (b ^^ e₆) (b ^^ e₇) (b ^^ e₈)) = b := by
  decide

theorem card_eq_sum_toNat (e : Finset (Fin 9)) : #e = ∑ i : Fin 9, (decide (i ∈ e)).toNat := by
  have h : ∀ i : Fin 9, (decide (i ∈ e)).toNat = if i ∈ e then 1 else 0 := fun i => by
    by_cases hi : i ∈ e <;> simp [hi]
  simp only [h, sum_boole, Nat.cast_id, filter_mem_eq_inter, univ_inter]

theorem received_eq_xor (b : Bool) (e : Finset (Fin 9)) (i : Fin 9) :
    received b e i = (b ^^ decide (i ∈ e)) := by
  by_cases hi : i ∈ e <;> simp [received, hi]

/-- Recursive majority on nine copies corrects every pattern of at most three flips. -/
theorem recursiveMajority_corrects (b : Bool) (e : Finset (Fin 9)) (he : #e ≤ 3) :
    recursiveMajority (received b e) = b := by
  rw [card_eq_sum_toNat] at he
  simp only [Fin.sum_univ_succ, Fin.sum_univ_zero] at he
  simp only [recursiveMajority, received_eq_xor]
  exact recursiveMajority_bits b _ _ _ _ _ _ _ _ _ (by simpa [add_assoc] using he)

/-- … and two flips in each of two groups defeat it with only four, which exact majority corrects. -/
theorem recursiveMajority_fails (b : Bool) :
    ∃ e : Finset (Fin 9), #e = 4 ∧ recursiveMajority (received b e) ≠ b :=
  ⟨{0, 1, 3, 4}, by decide, by cases b <;> decide⟩

end Pctc
