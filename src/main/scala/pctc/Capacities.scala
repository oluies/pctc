package pctc

import scotty.quantum.density.{Channel, ComplexMatrix, DensityMatrix}

/** One numerically estimated information measure. */
final case class InformationEstimate(bits: Double, point: Option[Vector[Double]], evaluations: Int):
  def isFinite: Boolean = bits.isFinite

/** The max-information and the Doeblin information of a channel, estimated by search.
  *
  * Both are taken from arXiv:2509.08965v3, where they give the asymptotic retrocausal capacity
  * `C_retro(N) = I_max(N) + I_doe^inf(N)` (Eq. (11)). `I_max` is additive, so its regularisation is
  * itself; `I_doe` is additive for the depolarizing channel (mictodiactic covariance,
  * Eq. (S196)-(S198)) and the erasure channel (Eq. (S215)), so for the channels swept here
  * `I_doe^inf = I_doe`. For the amplitude-damping channel `I_doe` is only a lower bound on its
  * regularisation, so `C_retro` as reported is an achievable rate rather than the capacity; the
  * paper makes the same point below Eq. (S226).
  *
  * Both definitions are semidefinite programs. Rather than pull in an SDP solver, the single-qubit
  * (or block-diagonal qutrit) operand is swept over a Bloch parametrisation and `lambda` is found by
  * bisection for each swept operator. The reported value is therefore an upper bound on the true
  * infimum, approached from above as the grid is refined.
  */
object Capacities:
  /** `N(Phi_{A'A})` on `(A', B)`: the Choi operator, normalised as a state. */
  def choiState(channel: Channel): DensityMatrix =
    require(channel.inputDims.length == 1, "The Choi state is built for a single-system channel.")

    channel.applyTo(DensityMatrix.maximallyEntangled(channel.inputDim), Seq(1))

  /** `pi_{A'} (x) operator`, with `pi` maximally mixed on a system of dimension `referenceDim`. */
  private def padded(referenceDim: Int, operator: ComplexMatrix): ComplexMatrix =
    ComplexMatrix.identity(referenceDim).scale(1d / referenceDim).tensor(operator)

  /** `I_max(N) = inf { log2 lambda : N(Phi_{A'A}) <= lambda pi_{A'} (x) sigma_B }`, Eq. (9).
    *
    * The infimum runs over real numbers and states of `B`.
    */
  def maxInformation(channel: Channel, grid: SearchGrid): InformationEstimate =
    val choi = choiState(channel).matrix
    val reference = channel.inputDim

    search(
      Bloch.stateFamily(channel.outputDim),
      grid,
      (operator, ceiling, tolerance) =>
        OperatorOrder.minLogLambda(choi, padded(reference, operator), ceiling, tolerance)
    )

  /** `I_doe(N) = inf { log2 lambda : pi_{A'} (x) tau_B <= lambda N(Phi_{A'A}) }`, Eq. (10).
    *
    * The infimum runs over real numbers and unit-trace Hermitian operators on `B`, which need not be
    * positive.
    */
  def doeblinInformation(channel: Channel, grid: SearchGrid): InformationEstimate =
    val choi = choiState(channel).matrix
    val reference = channel.inputDim

    search(
      Bloch.hermitianFamily(channel.outputDim, grid.radius),
      grid,
      (operator, ceiling, tolerance) =>
        OperatorOrder.minLogLambda(padded(reference, operator), choi, ceiling, tolerance)
    )

  /** `C_retro = I_max + I_doe`, infinite as soon as either part is. */
  def retrocausalCapacity(maxInformation: Double, doeblinInformation: Double): Double =
    maxInformation + doeblinInformation

  /** Locates the optimum with a tolerance that tightens as the search closes in, then evaluates the
    * winning point once at full precision. Scanning at full precision throughout would spend about
    * forty semidefinite tests on every grid point, almost all of them on points that lose.
    */
  private def search(
      family: OperatorFamily,
      grid: SearchGrid,
      objective: (ComplexMatrix, Double, Double) => Option[Double]
  ): InformationEstimate =
    val result = CoordinateSearch.minimize(
      family.box,
      grid,
      (point, ceiling, tolerance) => objective(family.build(point), ceiling, tolerance)
    )

    val refined = result.point.flatMap { point =>
      objective(family.build(point), OperatorOrder.Ceiling, OperatorOrder.LogTolerance)
    }

    val bits = refined.filter(_ < OperatorOrder.Ceiling).getOrElse(Double.PositiveInfinity)

    InformationEstimate(bits, result.point, result.evaluations + 1)
