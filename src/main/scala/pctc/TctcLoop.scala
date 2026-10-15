package pctc

import scotty.quantum.density.{ComplexMatrix, DensityMatrix}

/** The output of a transition-probability closed timelike curve, with its two terms kept apart.
  *
  * `coherent` is the post-selected contribution `P rho P^dagger` and `incoherent` is an ordinary
  * quantum channel; `weight` is the trace of their sum, which only normalises. Unlike a post-selected
  * curve there is no post-selection probability here, because a T-CTC never forbids a history.
  */
final case class TctcOutcome(
    normalized: DensityMatrix,
    coherent: DensityMatrix,
    incoherent: DensityMatrix,
    weight: Double
):
  /** How much of the output comes from the post-selected term rather than the ordinary channel. */
  def coherentFraction: Double = coherent.realTrace / weight

/** Allen's transition-probability closed timelike curve (T-CTC).
  *
  * J.-M. A. Allen, "Treating time travel quantum mechanically", Phys. Rev. A 90, 042107 (2014),
  * arXiv:1401.4933, defines the evolution as an average over chronology-violating pure states
  * weighted by their transition probability, Eq. (24):
  *
  * {{{
  *   rho_f = Z^-1 integral d[phi] U_phi rho_i U_phi^dagger,   U_phi = <phi|U|phi>
  * }}}
  *
  * Carrying out that integral over the Fubini-Study measure collapses it to Eq. (33), which is what
  * this implements:
  *
  * {{{
  *   rho_f = z [ P rho_i P^dagger + Tr_CV( U (rho_i (x) 1) U^dagger ) ],   P = Tr_CV(U)
  * }}}
  *
  * The first term is exactly the numerator of the post-selected map of [[PctcLoop]]; the second is an
  * ordinary, completely positive channel. Allen reads Eq. (33) as "a weighted mixture of the
  * corresponding P-CTC equation of motion with an ordinary quantum channel".
  *
  * That mixture is why T-CTCs behave differently from both earlier models. The second term is a
  * channel, so its trace cannot vanish for a non-zero input: **a T-CTC never assigns a history zero
  * probability**, where a post-selected curve forbids the grandfather paradox outright. And because
  * the two terms add incoherently, the output is generally mixed even when the post-selected one is
  * pure, which is what costs T-CTCs the ability to distinguish non-orthogonal states with certainty.
  */
object TctcLoop:
  /** Evaluates Eq. (33) for a unitary on the chronology-respecting system followed by the curve. */
  def evaluate(unitary: ComplexMatrix, crState: DensityMatrix, ctcDims: Vector[Int]): TctcOutcome =
    val ctcPositions = crState.dims.length until (crState.dims.length + ctcDims.length)

    require(
      unitary.rows == crState.dimension * ctcDims.product,
      s"U is ${unitary.rows}x${unitary.cols} but the joint system has dimension ${crState.dimension * ctcDims.product}"
    )

    // P rho P^dagger, the post-selected term.
    val p = PctcLoop.loopOperator(unitary, crState.dims ++ ctcDims, crState.dims.length)
    val coherent = DensityMatrix(crState.dims, p.multiply(crState.matrix).multiply(p.adjoint))

    // Tr_CV[ U (rho (x) 1) U^dagger ], the ordinary channel. The identity is deliberately not
    // normalised: Eq. (33) pairs `d Tr_CV[U (rho (x) 1/d) U^dagger]` with the term above.
    val padded = crState.tensor(DensityMatrix.identity(ctcDims))
    val evolved = DensityMatrix(padded.dims, unitary.multiply(padded.matrix).multiply(unitary.adjoint))
    val incoherent = evolved.partialTrace(ctcPositions)

    val total = coherent.plus(incoherent)
    val weight = total.realTrace

    require(weight > 0d, "A T-CTC always has a consistent solution, so this cannot vanish.")

    TctcOutcome(total.scale(1 / weight), coherent, incoherent, weight)

  /** The post-selected term on its own, for comparison with [[PctcLoop]]. */
  def postselectedTerm(unitary: ComplexMatrix, crState: DensityMatrix, ctcDims: Vector[Int]): DensityMatrix =
    evaluate(unitary, crState, ctcDims).coherent
