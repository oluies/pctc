package pctc

import scotty.quantum.density.{Channel, ComplexMatrix, DensityMatrix, Measurement}
import scotty.quantum.gate.StandardGate.RY
import scotty.quantum.math.Complex

/** The three demonstrations of the P-CTC framework, as pure computations.
  *
  * Nothing here prints or stores anything; [[Main]] does both.
  */
object Demos:
  private val trivial = DensityMatrix(Vector(1), ComplexMatrix.identity(1))
  private val qubitDim = 2

  /** `|psi>` for the message qubit used in the teleportation demonstration. */
  private val message = DensityMatrix.pure(Vector(2), Vector(Complex(0.6), Complex(0, 0.8)))

  /** Swaps `(E, B)` to `(F, A)`: the future end of the pair carries the message onward and the CTC
    * entrance receives it. This is the bipartite map whose loop is post-selected teleportation.
    */
  val swap: ComplexMatrix = ComplexMatrix.fromRows(
    Seq(
      Seq(Complex(1), Complex(0), Complex(0), Complex(0)),
      Seq(Complex(0), Complex(0), Complex(1), Complex(0)),
      Seq(Complex(0), Complex(1), Complex(0), Complex(0)),
      Seq(Complex(0), Complex(0), Complex(0), Complex(1))
    )
  )

  // -- (a) The grandfather paradox -------------------------------------------------------------

  /** Rotates the CTC qubit by `RY(theta)` inside an otherwise empty loop.
    *
    * The loop operator is `C = Tr[RY(theta)] = 2 cos(theta / 2)`, so the post-selection probability
    * is `|Tr RY(theta) / d|^2 = cos^2(theta / 2)`. At `theta = pi` the rotation is a bit flip, the
    * chronology-violating qubit would have to come back as its own negation, and the branch has
    * probability zero: the paradox is resolved by never happening.
    */
  def grandfather(steps: Int): Vector[GrandfatherRow] =
    require(steps >= 1, "The grandfather sweep needs at least one step.")

    Vector.tabulate(steps + 1): i =>
      val theta = Math.PI * i / steps
      val rotation = ComplexMatrix.fromRaw(RY(theta, 0).matrix)
      val t = PctcLoop.unitaryMap(rotation, Vector(1, qubitDim), Vector(1, qubitDim))
      val outcome = PctcLoop.evaluate(trivial.tensor(trivial), t, Channel.identity(qubitDim))

      GrandfatherRow(theta, outcome.postSelectionProbability, Math.pow(Math.cos(theta / 2), 2))

  // -- (b) Post-selected teleportation ---------------------------------------------------------

  /** Post-selected teleportation, first by hand and then as a P-CTC loop.
    *
    * By hand: the message is Bell-measured against the future end of a maximally entangled pair and
    * the `Phi+` outcome is post-selected. That outcome has probability 1/4, and conditioned on it the
    * past end carries the message exactly, with no correction to apply. Running the same bipartite
    * map inside the loop of Eq. (1) gives the same 1/4 as the denominator, which the loop then
    * divides out: the renormalised probability of the desirable outcome is 1.
    */
  def postselectedTeleportation: Vector[TeleportationRow] =
    val pair = DensityMatrix.maximallyEntangled(qubitDim)
    val t = PctcLoop.unitaryMap(swap, Vector(qubitDim, qubitDim), Vector(qubitDim, qubitDim))
    val loop = Channel.identity(qubitDim)

    // By hand, with a product message: subsystems are (M, past, future).
    val byHand = measureAndSteer(message.tensor(pair), Seq(0, 2), Seq(1))
    // By hand, with half of an entangled message: subsystems are (R, M, past, future).
    val entangled = measureAndSteer(pair.tensor(pair), Seq(1, 3), Seq(0, 2))

    Vector(
      TeleportationRow(
        "Bell measurement by hand, product message",
        byHand._1,
        byHand._1,
        overlapFidelity(message, byHand._2)
      ),
      TeleportationRow(
        "Bell measurement by hand, entangled message",
        entangled._1,
        entangled._1,
        overlapFidelity(pair, entangled._2)
      ),
      TeleportationRow(
        "P-CTC loop of Eq. (1), product message",
        loopProbability(trivial.tensor(message), t, loop),
        1d,
        overlapFidelity(message, PctcLoop.evaluate(trivial.tensor(message), t, loop).normalizedOrThrow.partialTrace(Seq(0)))
      ),
      TeleportationRow(
        "P-CTC loop of Eq. (1), entangled message",
        loopProbability(pair, t, loop),
        1d,
        overlapFidelity(pair, PctcLoop.evaluate(pair, t, loop).normalizedOrThrow)
      )
    )

  private def loopProbability(rhoRE: DensityMatrix, t: Channel, n: Channel): Double =
    PctcLoop.evaluate(rhoRE, t, n).postSelectionProbability

  /** Bell-measures `measured`, post-selects `Phi+`, and returns the probability together with the
    * normalised state left on `kept`.
    */
  private def measureAndSteer(
      state: DensityMatrix,
      measured: Seq[Int],
      kept: Seq[Int]
  ): (Double, DensityMatrix) =
    val outcome = Measurement.projectOnPureState(state, measured, Channel.bellBasis(0))

    (outcome.probability, outcome.normalized.partialTrace(measured))

  /** `Tr[target actual]`, which is the fidelity whenever `target` is pure, as it is here. */
  private def overlapFidelity(target: DensityMatrix, actual: DensityMatrix): Double =
    target.multiply(actual).realTrace

  // -- (c) Retrocausal messaging through a noisy CTC --------------------------------------------

  /** The max-information, the Doeblin information and their sum for one channel parameter.
    *
    * `postSelectionProbability` is the largest post-selection probability a
    * post-selected-teleportation loop through this noisy CTC can reach. By Lemma 1 / Eq. (S6) the
    * maximum achievable singlet fraction is `2^I_max / d_A^2`, and the loop of Eq. (1) contributes a
    * further factor `1 / d_A^2`, so the probability is `2^I_max / d_A^4`. For a noiseless CTC that is
    * `4 / 16 = 1/4`, which is exactly what demonstration (b) measures.
    */
  def retrocausal(spec: ChannelSpec, parameter: Double, grid: SearchGrid): RetrocausalRow =
    val channel = spec.build(parameter)
    val maxInformation = Capacities.maxInformation(channel, grid)
    val doeblin = Capacities.doeblinInformation(channel, grid)
    val dimension = channel.inputDim

    RetrocausalRow(
      channel = spec.name,
      parameterName = spec.parameterName,
      parameter = parameter,
      maxInformation = maxInformation.bits,
      doeblinInformation = doeblin.bits,
      retrocausalCapacity = Capacities.retrocausalCapacity(maxInformation.bits, doeblin.bits),
      maxInformationReference = spec.referenceMaxInformation(parameter),
      doeblinInformationReference = spec.referenceDoeblinInformation(parameter),
      retrocausalCapacityReference = spec.referenceRetrocausalCapacity(parameter),
      ordinaryClassicalCapacity = spec.ordinaryClassicalCapacity(parameter),
      postSelectionProbability = Math.pow(2d, maxInformation.bits) / Math.pow(dimension.toDouble, 4),
      evaluations = maxInformation.evaluations + doeblin.evaluations
    )

  /** The whole sweep: every channel family at every parameter. */
  def retrocausalSweep(parameters: Vector[Double], grid: SearchGrid): Vector[RetrocausalRow] =
    for spec <- ChannelSpec.all; parameter <- parameters yield retrocausal(spec, parameter, grid)
