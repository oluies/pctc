package pctc

import scotty.quantum.density.{ComplexMatrix, DensityMatrix, Eigen, NullSpace}
import scotty.quantum.math.Complex

/** The Deutsch-consistent states of a closed timelike curve: an affine family of density matrices.
  *
  * `directions` are traceless Hermitian operators, so every `at(coords)` has unit trace. The family
  * is a single point exactly when `directions` is empty.
  */
final case class ConsistentFamily(particular: DensityMatrix, directions: Vector[ComplexMatrix]):
  def isUnique: Boolean = directions.isEmpty

  def parameters: Int = directions.length

  def at(coordinates: Vector[Double]): DensityMatrix =
    require(coordinates.length == directions.length, "One coordinate per direction.")

    val shifted = coordinates.lazyZip(directions).foldLeft(particular.matrix): (acc, pair) =>
      acc.plus(pair._2.scale(pair._1))

    DensityMatrix(particular.dims, shifted)

/** What Deutsch's prescription predicts for one circuit. */
final case class DeutschOutcome(
    ctcState: DensityMatrix,
    output: DensityMatrix,
    family: ConsistentFamily,
    entropy: Double
):
  def isUnique: Boolean = family.isUnique

/** Deutsch's model of a closed timelike curve (D. Deutsch, Phys. Rev. D 44, 3197 (1991)).
  *
  * The chronology-violating system must come back as it left, which is the self-consistency
  * condition written as Eq. (1) of arXiv:1401.4933 (J.-M. A. Allen, Phys. Rev. A 90, 042107):
  *
  * {{{
  *   tau = D(tau) = Tr_CR[ U (rho_i (x) tau) U^dagger ]
  * }}}
  *
  * A solution always exists but need not be unique, and Deutsch resolves that with an
  * evolutionary principle selecting the solution of greatest von Neumann entropy.
  *
  * This is solved as a kernel, not by iteration. Iterating `D` does not generally converge: for the
  * grandfather paradox the map is `tau -> X tau X`, whose orbit has period two, so a fixed-point
  * loop oscillates forever and whatever it returns after a fixed number of steps is not a solution
  * at all. `(D - I)tau = 0` is exact, and the dimension of its kernel reports the non-uniqueness
  * directly instead of hiding it.
  *
  * Contrast [[PctcLoop]]: Deutsch's condition is always satisfiable, at the cost of mixing the
  * chronology-violating system, so nothing is ever forbidden. A post-selected curve instead gives
  * the paradoxical branch probability zero.
  */
object DeutschCtc:
  /** `tau -> Tr_CR[ U (rho_i (x) tau) U^dagger ]`, the map whose fixed points are the solutions. */
  def consistencyMap(unitary: ComplexMatrix, crState: DensityMatrix, ctcDims: Vector[Int])(
      ctcState: DensityMatrix
  ): DensityMatrix =
    val full = crState.tensor(ctcState)

    require(unitary.rows == full.dimension, s"U is ${unitary.rows}x${unitary.cols}, state is ${full.dimension}")

    val evolved = DensityMatrix(full.dims, unitary.multiply(full.matrix).multiply(unitary.adjoint))

    evolved.partialTrace(crState.dims.indices)

  /** The chronology-respecting output, `Tr_CTC[ U (rho_i (x) tau) U^dagger ]`. */
  def output(unitary: ComplexMatrix, crState: DensityMatrix, ctcState: DensityMatrix): DensityMatrix =
    val full = crState.tensor(ctcState)
    val evolved = DensityMatrix(full.dims, unitary.multiply(full.matrix).multiply(unitary.adjoint))
    val ctcPositions = crState.dims.length until full.dims.length

    evolved.partialTrace(ctcPositions)

  /** The consistency map as a `d^2 x d^2` matrix acting on row-major stacked operators. */
  def superoperator(unitary: ComplexMatrix, crState: DensityMatrix, ctcDims: Vector[Int]): ComplexMatrix =
    val d = ctcDims.product
    val map = consistencyMap(unitary, crState, ctcDims)

    // Column (i, j) is the map applied to the basis operator |i><j|, stacked the same way.
    val columns = for i <- 0 until d; j <- 0 until d yield
      val unit = ComplexMatrix(d, d, Vector.tabulate(d * d)(k => if k == i * d + j then Complex(1) else Complex(0)))

      map(DensityMatrix(ctcDims, unit)).matrix.data

    ComplexMatrix(d * d, d * d, Vector.tabulate(d * d * d * d)(k => columns(k % (d * d))(k / (d * d))))

  /** Every Deutsch-consistent CTC state, as an affine family. */
  def consistentStates(
      unitary: ComplexMatrix,
      crState: DensityMatrix,
      ctcDims: Vector[Int],
      tolerance: Double = 1e-9
  ): ConsistentFamily =
    val d = ctcDims.product
    val kernel = NullSpace.basis(superoperator(unitary, crState, ctcDims).minus(ComplexMatrix.identity(d * d)), tolerance)

    require(kernel.nonEmpty, "The consistency condition has no solution, which cannot happen for a channel.")

    // The map preserves adjoints, so the Hermitian and anti-Hermitian parts of a fixed operator are
    // themselves fixed. Taking both gives a real-linear basis of the fixed Hermitian operators.
    val hermitian = kernel.flatMap: vector =>
      val m = ComplexMatrix(d, d, vector)

      Vector(m.plus(m.adjoint).scale(0.5), m.minus(m.adjoint).scale(Complex(0, -0.5)))

    val independent = orthogonalize(hermitian.filter(_.data.exists(c => Complex.abs(c) > tolerance)), tolerance)
    val traced = independent.find(m => Math.abs(m.trace.r) > tolerance)

    require(traced.nonEmpty, "No unit-trace solution, which cannot happen for a trace-preserving map.")

    val particular = traced.get.scale(1 / traced.get.trace.r)
    val directions = orthogonalize(
      independent.filterNot(_ eq traced.get).map(m => m.minus(particular.scale(m.trace.r))),
      tolerance
    ).map(m => m.scale(1 / frobenius(m)))

    ConsistentFamily(DensityMatrix(ctcDims, particular), directions)

  /** Deutsch's maximum-entropy solution, and the chronology-respecting output it induces. */
  def solve(
      unitary: ComplexMatrix,
      crState: DensityMatrix,
      ctcDims: Vector[Int],
      grid: SearchGrid = SearchGrid(steps = 17, rounds = 40, radius = 1d)
  ): DeutschOutcome =
    val family = consistentStates(unitary, crState, ctcDims)

    val chosen =
      if family.isUnique then family.particular
      else
        val box = Vector.fill(family.parameters)((-2d, 2d))
        val result = CoordinateSearch.minimize(
          box,
          grid,
          (point, ceiling, _) =>
            val candidate = family.at(point)

            Option.when(candidate.isPositiveSemidefinite(1e-9)) {
              -entropy(candidate)
            }.filter(_ < ceiling)
        )

        result.point.map(family.at).getOrElse(family.particular)

    DeutschOutcome(chosen, output(unitary, crState, chosen), family, entropy(chosen))

  /** Von Neumann entropy in bits. */
  def entropy(state: DensityMatrix): Double =
    state.eigenvalues.foldLeft(0d): (sum, lambda) =>
      if lambda <= 1e-12 then sum else sum - lambda * Math.log(lambda) / Math.log(2)

  private def frobenius(m: ComplexMatrix): Double =
    Math.sqrt(m.data.foldLeft(0d)((s, c) => s + c.r * c.r + c.i * c.i))

  /** Gram-Schmidt under the Hilbert-Schmidt inner product, dropping dependent elements. */
  private def orthogonalize(ms: Vector[ComplexMatrix], tolerance: Double): Vector[ComplexMatrix] =
    ms.foldLeft(Vector.empty[ComplexMatrix]): (kept, candidate) =>
      val reduced = kept.foldLeft(candidate): (acc, basis) =>
        val overlap = basis.adjoint.multiply(acc).trace.r / (frobenius(basis) * frobenius(basis))

        acc.minus(basis.scale(overlap))

      if frobenius(reduced) > tolerance * Math.max(1d, frobenius(candidate)) then kept :+ reduced else kept
