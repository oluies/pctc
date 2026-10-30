package pctc

import scotty.quantum.density.{Channel, ComplexMatrix, DensityMatrix}
import scotty.quantum.gate.StandardGate.{H, X}
import scotty.quantum.math.Complex

/** One row of the repetition-code side demonstration. `simulated` comes from the density matrix,
  * `closedForm` from the formula it has to match.
  */
final case class RepetitionRow(noise: String, distance: Int, probability: Double, simulated: Double, closedForm: Double)

/** A side demonstration, not part of arXiv:2509.08965: why cat-qubit designs leave out protection
  * against bit flips.
  *
  * Alice & Bob's estimate of 22 593 cat qubits for RSA-2048 (IACR ePrint 2026/2665) assumes qubits
  * whose bit flips are rare enough to ignore, a bit-flip time of about a hundred years, so that a
  * repetition code against phase flips alone is enough. The paper's own Appendix C9 puts today's
  * measured cat qubits more than eight orders of magnitude short of that. The phase-flip code of distance `d` stores logical `|0>` as `|+>^d` and logical `|1>` as
  * `|->^d`, and reads a logical value by majority vote over `X`-basis measurements.
  *
  * This is the code-capacity model: noise once on every qubit, then perfect measurement. The paper's
  * figures come from circuit-level noise during repeated syndrome extraction, which is a different
  * and harder model; nothing here reproduces them.
  */
object RepetitionCode:
  private val Sqrt2Inverse = 1 / Math.sqrt(2)

  private val plus = DensityMatrix.pure(Vector(2), Seq(Complex(Sqrt2Inverse), Complex(Sqrt2Inverse)))
  private val hadamard = Channel.fromMatrix(H(0).matrix)

  /** `rho -> (1 - p) rho + p X rho X`. */
  def bitFlip(p: Double): Channel =
    require(p >= 0d && p <= 1d, "A bit-flip probability has to lie in [0, 1].")

    Channel(
      Vector(2),
      Vector(2),
      Vector(
        ComplexMatrix.identity(2).scale(Math.sqrt(1 - p)),
        ComplexMatrix.fromRaw(X(0).matrix).scale(Math.sqrt(p))
      )
    )

  private def requireDistance(distance: Int): Unit =
    require(distance >= 1 && distance % 2 == 1, "A repetition code needs an odd distance.")

  private def onEveryQubit(rho: DensityMatrix, channel: Channel): DensityMatrix =
    rho.dims.indices.foldLeft(rho)((state, qubit) => channel.applyTo(state, Seq(qubit)))

  /** Logical `|0>` dephases with probability `p` on every qubit and is read by majority vote.
    *
    * After a Hadamard on every qubit, basis state `i` is the measurement outcome with `|->` wherever
    * `i` has a one, so the decoder fails exactly on the outcomes with more than `d / 2` ones.
    */
  def phaseFlipFailure(distance: Int, p: Double): Double =
    requireDistance(distance)

    val encoded = Vector.fill(distance - 1)(plus).foldLeft(plus)(_.tensor(_))
    val measured = onEveryQubit(onEveryQubit(encoded, Channel.dephasing(p)), hadamard)

    (0 until measured.dimension).filter(i => Integer.bitCount(i) > distance / 2).map(i => measured.matrix.realAt(i, i)).sum

  /** `sum_{k > d/2} C(d, k) p^k (1 - p)^(d - k)`: more than half the qubits flipped. Proved, as the
    * exact failure probability of majority vote, by `Pctc.failureProbability_eq` in
    * `lean/Pctc/RepetitionCode.lean`.
    */
  def phaseFlipFailureClosedForm(distance: Int, p: Double): Double =
    requireDistance(distance)

    (distance / 2 + 1 to distance).map(k => binomial(distance, k) * Math.pow(p, k) * Math.pow(1 - p, distance - k)).sum

  /** Logical `|+> = (|+>^d + |->^d) / sqrt 2` suffers a bit flip with probability `p` on every qubit.
    *
    * An `X` leaves `|+>` alone and multiplies `|->` by -1, so on this code one bit flip is already a
    * logical error: it takes logical `|+>` to logical `|->`. The `X`-basis checks commute with it and
    * cannot see it, so the probability returned, of ending in logical `|->`, is the decoder's failure
    * rate whatever it does.
    */
  def bitFlipFailure(distance: Int, p: Double): Double =
    requireDistance(distance)

    val noisy = onEveryQubit(DensityMatrix.pure(Vector.fill(distance)(2), logical(distance, 1d)), bitFlip(p))
    val minus = ComplexMatrix.column(logical(distance, -1d))

    minus.adjoint.multiply(noisy.matrix).multiply(minus).realAt(0, 0)

  /** `(1 - (1 - 2p)^d) / 2`, an odd number of flips; about `d p` when `p` is small, so it grows with
    * the distance. Written with `expm1` and `log1p`, which keep it accurate for small `p`.
    */
  def bitFlipFailureClosedForm(distance: Int, p: Double): Double =
    requireDistance(distance)

    if p >= 0.5 then (1 - Math.pow(1 - 2 * p, distance)) / 2 else -Math.expm1(distance * Math.log1p(-2 * p)) / 2

  /** The amplitudes of `(|+>^d + sign |->^d) / sqrt 2`, given that `|+>^d` has amplitude `2^(-d/2)`
    * everywhere and `|->^d` has `(-1)^popcount(i) 2^(-d/2)`.
    */
  private def logical(distance: Int, sign: Double): Vector[Complex] =
    val scale = Sqrt2Inverse / Math.pow(2, distance / 2d)

    Vector.tabulate(1 << distance): i =>
      Complex(scale * (1 + sign * (if Integer.bitCount(i) % 2 == 0 then 1 else -1)))

  private def binomial(n: Int, k: Int): Double =
    (1 to k).foldLeft(1d)((acc, i) => acc * (n - k + i) / i)

  /** The paper's phase-flip probability per cycle. */
  val PaperPhaseFlip: Double = 5e-4

  /** Both noise models at distances 1 to 7 and three probabilities, including the paper's `p_Z`. */
  def rows: Vector[RepetitionRow] =
    val distances = Vector(1, 3, 5, 7)
    val probabilities = Vector(0.1, 0.01, PaperPhaseFlip)

    val phase = for p <- probabilities; d <- distances yield
      RepetitionRow("phase flip", d, p, phaseFlipFailure(d, p), phaseFlipFailureClosedForm(d, p))
    val bit = for p <- probabilities; d <- distances yield
      RepetitionRow("bit flip", d, p, bitFlipFailure(d, p), bitFlipFailureClosedForm(d, p))

    phase ++ bit
