package pctc

import scotty.quantum.density.{Channel, ComplexMatrix, DensityMatrix}

/** One evaluation of the P-CTC loop map.
  *
  * `unnormalized` is the numerator of Eq. (1), `postSelectionProbability` its denominator. Both are
  * needed: the probability is what a loop renormalises away, and when it vanishes (the grandfather
  * paradox at its extreme) there is no normalised output at all.
  */
final case class LoopOutcome(unnormalized: DensityMatrix, postSelectionProbability: Double):
  def normalized(threshold: Double = 1e-12): Option[DensityMatrix] =
    Option.when(postSelectionProbability > threshold)(unnormalized.normalized)

  def normalizedOrThrow: DensityMatrix = normalized().getOrElse(
    throw new IllegalStateException(
      s"The post-selected branch has probability $postSelectionProbability; there is no normalised output."
    )
  )

/** The nonlinear map a noisy post-selected closed timelike curve induces on a bipartite channel.
  *
  * This is Eq. (1) of arXiv:2509.08965v3 (Eq. (2) of the first arXiv version), "Retrocausal capacity
  * of a quantum channel: Communicating through noisy closed timelike curves", Ji, Lloyd and Wilde:
  *
  * {{{
  *                Tr_{AA'}[ Phi_{AA'} T_{EB->FA}[ rho_RE (x) N_{A->B}[Phi_{AA'}] ] ]
  *   rho_RE  |->  -------------------------------------------------------------------
  *                  Tr[ Phi_{AA'} T_{EB->FA}[ rho_RE (x) N_{A->B}[Phi_{AA'}] ] ]
  * }}}
  *
  * with `R` an arbitrary reference system and `Phi_{AA'}` the standard maximally entangled state of
  * `A` and `A'`. `N_{A->B}` is the noisy P-CTC: `B` precedes `A` in time. `T_{EB->FA}` is any
  * completely positive map, which is why [[scotty.quantum.density.Channel]] does not insist on trace
  * preservation.
  *
  * Nothing here sends information backward in time. This is a numerical model of the paper's
  * framework and nothing more.
  */
object PctcLoop:
  /** Evaluates the map above.
    *
    * @param rhoRE two subsystems, the reference `R` first and the map's input `E` second. Use a
    *              dimension of 1 for either when it is not needed.
    * @param t     the bipartite map, with input `(E, B)` and output `(F, A)` in that order.
    * @param n     the noisy P-CTC channel `A -> B`.
    */
  def evaluate(rhoRE: DensityMatrix, t: Channel, n: Channel): LoopOutcome =
    require(rhoRE.dims.length == 2, "rho_RE has to factorise as R (x) E.")
    require(t.inputDims.length == 2, "T has to have input factors (E, B).")
    require(t.outputDims.length == 2, "T has to have output factors (F, A).")

    val Vector(dR, dE) = rhoRE.dims
    val Vector(dF, dA) = t.outputDims

    require(t.inputDims == Vector(dE, n.outputDim), s"T's input ${t.inputDims} has to be (E, B).")
    require(dA == n.inputDim, s"T's output system A has dimension $dA but N consumes ${n.inputDim}.")

    val phi = DensityMatrix.maximallyEntangled(dA)         // on (A, A')
    val loopResource = n.applyTo(phi, Seq(0))              // on (B, A')
    val input = rhoRE.tensor(loopResource)                 // on (R, E, B, A')
    val evolved = t.applyTo(input, Seq(1, 2))              // on (R, F, A, A')

    val projector = ComplexMatrix.identity(dR * dF).tensor(phi.matrix)
    val weighted = DensityMatrix(evolved.dims, projector.multiply(evolved.matrix))
    val numerator = weighted.partialTrace(Seq(2, 3))       // on (R, F)

    LoopOutcome(numerator, numerator.realTrace)

  /** `T(X) = U X U^dagger` for a `U` taking `(E, B)` to `(F, A)`. */
  def unitaryMap(unitary: ComplexMatrix, input: Vector[Int], output: Vector[Int]): Channel =
    Channel(input, output, Vector(unitary))

  /** The Lloyd et al. loop operator `C = Tr_CTC[U]` (arXiv:1005.2219).
    *
    * For a unitary `U` on the chronology-respecting system tensored with the chronology-violating
    * (CTC) system, the noiseless P-CTC acts on the chronology-respecting system as
    * `rho -> C rho C^dagger / Tr[C rho C^dagger]`, where `C` is `U` partially traced over the CTC
    * factor. The loop map above reproduces exactly that, with a post-selection probability of
    * `Tr[C rho C^dagger] / d_CTC^2`.
    */
  def loopOperator(unitary: ComplexMatrix, dims: Vector[Int], ctcSubsystem: Int): ComplexMatrix =
    require(unitary.isSquare, "The loop unitary has to be square.")

    DensityMatrix(dims, unitary).partialTrace(Seq(ctcSubsystem)).matrix
