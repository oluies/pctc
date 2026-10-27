/-
Reproducible check that the results carry no `sorry`.

    cd lean && lake env lean check_axioms.lean

Each theorem should report at most [propext, Classical.choice, Quot.sound], the three
standard axioms of Lean's logic. The appearance of `sorryAx` anywhere would mean a gap.
-/
import Pctc.Loop
import Pctc.BlockDiagonal

#print axioms Pctc.loopInput_apply
#print axioms Pctc.loopEvolved_apply
#print axioms Pctc.loopNumerator_eq
#print axioms Pctc.trace_loopNumerator
#print axioms Pctc.loop_normalized_eq
#print axioms Pctc.grandfather_postSelection
#print axioms Pctc.grandfather_postSelection_pi
#print axioms Pctc.maxInformation_pinch
#print axioms Pctc.doeblinInformation_pinch
#print axioms Pctc.choi_invariant
#print axioms Pctc.erasureChoi_kept
#print axioms Pctc.erasureChoi_flag
#print axioms Pctc.erasureChoi_invariant
#print axioms Pctc.erasure_maxInformation_blockDiagonal
#print axioms Pctc.erasure_doeblinInformation_blockDiagonal
