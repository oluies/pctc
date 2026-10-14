package pctc

import scotty.quantum.density.{ComplexMatrix, Eigen}
import scotty.quantum.math.Complex

/** The semidefinite order tests behind the max- and Doeblin informations.
  *
  * Both quantities of arXiv:2509.08965 are infima of `log2 lambda` subject to an operator
  * inequality, so everything reduces to deciding `X <= 2^l Y` for Hermitian `X` and positive
  * semidefinite `Y`, and then finding the smallest `l` for which that holds.
  *
  * Because `Y >= 0`, feasibility is monotone in `l`: if `l` works then so does every larger value.
  * That makes plain bisection on `log2 lambda` both correct and robust in the presence of a singular
  * `Y`, which an explicit generalised-eigenvalue formula would have to special-case.
  */
object OperatorOrder:
  /** Smallest `log2 lambda` considered. Below this the two operators are treated as equal. */
  val Floor: Double = -40d

  /** Largest `log2 lambda` considered. Above this the infimum is reported as infinite. */
  val Ceiling: Double = 32d

  /** Bisection tolerance on `log2 lambda`, in bits. */
  val LogTolerance: Double = 1e-10

  /** Relative tolerance for the eigenvalue sign test, which absorbs Jacobi round-off.
    *
    * This has to stay well below `2 ^ -Ceiling`. When no `lambda` at all satisfies the inequality,
    * the violation on the kernel of `reference` shrinks like `2 ^ -logLambda` under the
    * normalisation used by [[dominates]], so a tolerance above `2 ^ -Ceiling` would accept the
    * ceiling and report a finite value where the honest answer is infinite.
    */
  val PositivityTolerance: Double = 1e-12

  /** Is `candidate <= 2^logLambda * reference`?
    *
    * For a non-negative exponent the inequality is divided by `lambda` rather than multiplied, so the
    * entries being compared never grow with `lambda` and a single relative tolerance stays
    * meaningful across the whole bisection range.
    */
  def dominates(candidate: ComplexMatrix, reference: ComplexMatrix, logLambda: Double): Boolean =
    val difference =
      if logLambda >= 0d then reference.minus(candidate.scale(Math.pow(2d, -logLambda)))
      else reference.scale(Math.pow(2d, logLambda)).minus(candidate)

    // The tolerance is relative to the matrix actually being tested. Taking it relative to the raw
    // operands instead would let a `tau_B` with a long Bloch vector inflate it until a genuinely
    // infeasible inequality passed.
    val scale = difference.data.foldLeft(1e-300)((m, c) => Math.max(m, Complex.abs(c)))

    Eigen.minEigenvalue(difference) >= -PositivityTolerance * scale

  /** `inf { log2 lambda : candidate <= lambda * reference }`, or `+inf` when no `lambda` works.
    *
    * `ceiling` doubles as an incumbent: pass the best value found so far and this returns `None`
    * immediately for any pair that cannot beat it, which is what keeps the grid searches cheap.
    */
  def minLogLambda(
      candidate: ComplexMatrix,
      reference: ComplexMatrix,
      ceiling: Double = Ceiling
  ): Option[Double] =
    if !dominates(candidate, reference, ceiling) then None
    else Some(bisect(candidate, reference, Floor, ceiling))

  private def bisect(
      candidate: ComplexMatrix,
      reference: ComplexMatrix,
      low: Double,
      high: Double
  ): Double =
    if high - low <= LogTolerance then high
    else
      val middle = low + (high - low) / 2

      if dominates(candidate, reference, middle) then bisect(candidate, reference, low, middle)
      else bisect(candidate, reference, middle, high)
