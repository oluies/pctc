package pctc

import scotty.quantum.density.{Cholesky, ComplexMatrix}

/** The semidefinite order tests behind the max- and Doeblin informations.
  *
  * Both quantities of arXiv:2509.08965 are infima of `log2 lambda` subject to an operator
  * inequality, so everything reduces to deciding `X <= 2^l Y` for Hermitian `X` and positive
  * semidefinite `Y`, and then finding the smallest `l` for which that holds.
  *
  * Because `Y >= 0`, feasibility is monotone in `l`: if `l` works then so does every larger value.
  * That makes plain bisection on `log2 lambda` both correct and robust in the presence of a singular
  * `Y`, which an explicit generalised-eigenvalue formula would have to special-case.
  *
  * The decision itself goes through [[Cholesky]] rather than [[scotty.quantum.density.Eigen]]: the
  * question is only whether the smallest eigenvalue clears a threshold, and a factorisation answers
  * that without computing any eigenvalue.
  */
object OperatorOrder:
  /** Smallest `log2 lambda` considered. Below this the two operators are treated as equal. */
  val Floor: Double = -40d

  /** Largest `log2 lambda` considered. Above this the infimum is reported as infinite. */
  val Ceiling: Double = 32d

  /** Bisection tolerance on `log2 lambda`, in bits, for a fully refined answer. */
  val LogTolerance: Double = 1e-10

  /** Relative tolerance for the positivity test, which absorbs round-off.
    *
    * This has to stay well below `2 ^ -Ceiling`. When no `lambda` at all satisfies the inequality,
    * the violation on the kernel of `reference` shrinks like `2 ^ -logLambda` under the
    * normalisation used by [[dominates]], so a tolerance above `2 ^ -Ceiling` would accept the
    * ceiling and report a finite value where the honest answer is infinite. It also has to stay
    * well above the factorisation's own round-off, which is around `1e-15` relative.
    */
  val PositivityTolerance: Double = 1e-12

  /** Is `candidate <= 2^logLambda * reference`?
    *
    * For a non-negative exponent the inequality is divided by `lambda` rather than multiplied, so the
    * entries being compared never grow with `lambda` and a single relative tolerance stays
    * meaningful across the whole bisection range.
    */
  def dominates(candidate: ComplexMatrix, reference: ComplexMatrix, logLambda: Double): Boolean =
    // The tolerance is relative to the difference actually being tested. Taking it relative to the
    // raw operands instead would let a `tau_B` with a long Bloch vector inflate it until a genuinely
    // infeasible inequality passed.
    if logLambda >= 0d then
      Cholesky.differenceIsPositiveDefinite(reference, candidate, Math.pow(2d, -logLambda), PositivityTolerance)
    else
      Cholesky.differenceIsPositiveDefinite(
        reference.scale(Math.pow(2d, logLambda)),
        candidate,
        1d,
        PositivityTolerance
      )

  /** `inf { log2 lambda : candidate <= lambda * reference }`, or `None` when no `lambda` works.
    *
    * `ceiling` doubles as an incumbent: pass the best value found so far and this returns `None`
    * after a single semidefinite test for any pair that cannot beat it, which is what keeps the grid
    * searches cheap. `tolerance` is how tightly the answer is bisected, in bits; a search scanning
    * for the argument of the minimum can afford a loose one and refine the winner afterwards.
    */
  def minLogLambda(
      candidate: ComplexMatrix,
      reference: ComplexMatrix,
      ceiling: Double = Ceiling,
      tolerance: Double = LogTolerance
  ): Option[Double] =
    if !dominates(candidate, reference, ceiling) then None
    else Some(bisect(candidate, reference, Floor, ceiling, tolerance))

  private def bisect(
      candidate: ComplexMatrix,
      reference: ComplexMatrix,
      low: Double,
      high: Double,
      tolerance: Double
  ): Double =
    if high - low <= tolerance then high
    else
      val middle = low + (high - low) / 2

      if dominates(candidate, reference, middle) then bisect(candidate, reference, low, middle, tolerance)
      else bisect(candidate, reference, middle, high, tolerance)
