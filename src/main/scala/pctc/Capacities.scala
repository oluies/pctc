package pctc

import scotty.quantum.density.{Channel, ComplexMatrix, DensityMatrix}

/** The two information measures of one channel, computed together. */
final case class ChannelInformations(maxInformation: InformationEstimate, doeblinInformation: InformationEstimate)

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

  /** `pi_{A'}`, maximally mixed on a system of the given dimension.
    *
    * Built once per search and closed over, not per swept operator: the searches evaluate their
    * objective a few million times, and rebuilding a constant there is a few million allocations
    * thrown away.
    */
  private def maximallyMixed(dimension: Int): ComplexMatrix =
    ComplexMatrix.identity(dimension).scale(1d / dimension)

  /** `I_max(N) = inf { log2 lambda : N(Phi_{A'A}) <= lambda pi_{A'} (x) sigma_B }`, Eq. (9).
    *
    * The infimum runs over real numbers and states of `B`.
    */
  def maxInformation(channel: Channel, grid: SearchGrid): InformationEstimate =
    maxInformationOf(choiState(channel).matrix, maximallyMixed(channel.inputDim), channel, grid)

  /** `I_doe(N) = inf { log2 lambda : pi_{A'} (x) tau_B <= lambda N(Phi_{A'A}) }`, Eq. (10).
    *
    * The infimum runs over real numbers and unit-trace Hermitian operators on `B`, which need not be
    * positive.
    */
  def doeblinInformation(channel: Channel, grid: SearchGrid): InformationEstimate =
    doeblinInformationOf(choiState(channel).matrix, maximallyMixed(channel.inputDim), channel, grid)

  /** Both informations, sharing one Choi operator and one `pi` between the two searches.
    *
    * Each search is still run exactly once. Asking for the two measures separately is also fine; it
    * only recomputes the Choi operator, which is cheap next to either search.
    */
  def informations(channel: Channel, grid: SearchGrid): ChannelInformations =
    val choi = choiState(channel).matrix
    val pi = maximallyMixed(channel.inputDim)

    ChannelInformations(
      maxInformation = maxInformationOf(choi, pi, channel, grid),
      doeblinInformation = doeblinInformationOf(choi, pi, channel, grid)
    )

  private def maxInformationOf(
      choi: ComplexMatrix,
      pi: ComplexMatrix,
      channel: Channel,
      grid: SearchGrid
  ): InformationEstimate =
    search(
      Bloch.stateFamily(channel.outputDim),
      grid,
      (operator, bound) => OperatorOrder.bracket(choi, pi.tensor(operator), bound)
    )

  private def doeblinInformationOf(
      choi: ComplexMatrix,
      pi: ComplexMatrix,
      channel: Channel,
      grid: SearchGrid
  ): InformationEstimate =
    search(
      Bloch.hermitianFamily(channel.outputDim, grid.radius),
      grid,
      (operator, bound) => OperatorOrder.bracket(pi.tensor(operator), choi, bound)
    )

  /** `C_retro = I_max + I_doe`, infinite as soon as either part is. */
  def retrocausalCapacity(maxInformation: Double, doeblinInformation: Double): Double =
    maxInformation + doeblinInformation

  /** Locates the optimum, then settles the winning bracket to full precision.
    *
    * The search itself refines only as far as its comparisons demand, so a point that loses costs a
    * single semidefinite test rather than a full bisection, and the winner arrives already partly
    * refined.
    */
  private def search(
      family: OperatorFamily,
      grid: SearchGrid,
      objective: (ComplexMatrix, Double) => Option[Bracket]
  ): InformationEstimate =
    val result = CoordinateSearch
      .minimize(family.box, grid, (point, bound) => objective(family.build(point), bound))
      .narrowed(OperatorOrder.LogTolerance)

    val bits =
      if result.value >= OperatorOrder.Ceiling then Double.PositiveInfinity else result.value

    InformationEstimate(bits, result.point, result.evaluations)
