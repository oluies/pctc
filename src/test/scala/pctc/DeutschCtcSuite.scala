package pctc

import scotty.quantum.density.{Channel, ComplexMatrix, DensityMatrix}
import scotty.quantum.math.Complex

/** Checked against the worked examples published in J.-M. A. Allen, "Treating time travel quantum
  * mechanically", Phys. Rev. A 90, 042107 (2014), arXiv:1401.4933. Those are solutions someone else
  * derived for circuits someone else chose, which is the point: the rest of this suite would pass
  * just as happily on a model I had got subtly wrong.
  */
class DeutschCtcSuite extends munit.FunSuite:
  private val tolerance = 1e-9
  private val qubit = Vector(2)

  private def ketBra(amplitudes: Seq[Complex]): DensityMatrix =
    DensityMatrix.pure(Vector.fill(Math.round(Math.log(amplitudes.length.toDouble) / Math.log(2)).toInt)(2), amplitudes)

  private val zero = DensityMatrix.basis(qubit, 0)
  private val one = DensityMatrix.basis(qubit, 1)
  private val half = DensityMatrix.maximallyMixed(qubit)
  private val inverseRoot2 = 1 / Math.sqrt(2)
  private val minus = ketBra(Seq(Complex(inverseRoot2), Complex(-inverseRoot2)))
  private val plus = ketBra(Seq(Complex(inverseRoot2), Complex(inverseRoot2)))

  // -- The grandfather paradox ------------------------------------------------------------------

  test("the grandfather CNOT has a non-unique solution whose maximum-entropy member is I/2") {
    val u = Circuits.cnot(qubits = 2, control = 0, target = 1)
    val outcome = DeutschCtc.solve(u, one, qubit)

    assert(!outcome.isUnique, "Deutsch's condition is famously degenerate here")
    assertEquals(outcome.family.parameters, 1)
    assert(outcome.ctcState.approximatelyEquals(half, 1e-6), outcome.ctcState.toHumanString)
    assertEqualsDouble(outcome.entropy, 1d, 1e-6)
    // The survival probability Deutsch predicts, which a naive fixed-point iteration misses.
    assertEqualsDouble(outcome.ctcState(1, 1).r, 0.5, 1e-6)
  }

  test("every member of the grandfather family really is a fixed point") {
    val u = Circuits.cnot(qubits = 2, control = 0, target = 1)
    val family = DeutschCtc.consistentStates(u, one, qubit)
    val map = DeutschCtc.consistencyMap(u, one, qubit)

    for coordinate <- Seq(-0.9, -0.3, 0d, 0.4, 1.0) do
      val candidate = family.at(Vector(coordinate))

      assert(map(candidate).approximatelyEquals(candidate, tolerance), s"coordinate $coordinate")
  }

  test("a loop with no interaction leaves the curve unconstrained, with the identity as its map") {
    val u = ComplexMatrix.identity(4)
    val family = DeutschCtc.consistentStates(u, zero, qubit)

    // Every state is consistent, so the family spans the whole Bloch ball: three directions.
    assertEquals(family.parameters, 3)
    assertEqualsDouble(DeutschCtc.solve(u, zero, qubit).entropy, 1d, 1e-6)
  }

  // -- Allen, the unproven-theorem circuit -----------------------------------------------------

  /** `U = SWAP_MT CNOT_BM CNOT_TB` on qubits (B, M, T), with T the chronology-violating one. */
  private val unprovenTheorem =
    Circuits.inOrder(
      Circuits.cnot(qubits = 3, control = 2, target = 0),
      Circuits.cnot(qubits = 3, control = 0, target = 1),
      Circuits.swap(qubits = 3, first = 1, second = 2)
    )

  private val bothZero = DensityMatrix.basis(Vector(2, 2), 0)

  test("Allen Eq. (4): the unproven-theorem circuit has a one-parameter family of solutions") {
    val family = DeutschCtc.consistentStates(unprovenTheorem, bothZero, qubit)

    assertEquals(family.parameters, 1)

    // The paper gives the family as alpha|0><0| + (1 - alpha)|1><1|.
    val map = DeutschCtc.consistencyMap(unprovenTheorem, bothZero, qubit)

    for alpha <- Seq(0d, 0.25, 0.5, 0.75, 1d) do
      val tau = DensityMatrix(qubit, zero.matrix.scale(alpha).plus(one.matrix.scale(1 - alpha)))

      assert(map(tau).approximatelyEquals(tau, tolerance), s"alpha = $alpha is in the published family")
  }

  test("Allen Eq. (4): the maximum-entropy rule picks alpha = 1/2 and a classically correlated output") {
    val outcome = DeutschCtc.solve(unprovenTheorem, bothZero, qubit)

    assert(outcome.ctcState.approximatelyEquals(half, 1e-6), outcome.ctcState.toHumanString)

    val expected = DensityMatrix(
      Vector(2, 2),
      DensityMatrix.basis(Vector(2, 2), 0).matrix.scale(0.5)
        .plus(DensityMatrix.basis(Vector(2, 2), 3).matrix.scale(0.5))
    )

    assert(outcome.output.approximatelyEquals(expected, 1e-6), outcome.output.toHumanString)
  }

  test("Allen Eq. (9): a post-selected curve on the same circuit gives a pure entangled output") {
    // P = Tr_CV(U), with the chronology-violating qubit last.
    val p = PctcLoop.loopOperator(unprovenTheorem, Vector(2, 2, 2), ctcSubsystem = 2)
    val input = ComplexMatrix.column(Vector(Complex(1), Complex(0), Complex(0), Complex(0)))
    val out = p.multiply(input)
    val weight = out.adjoint.multiply(out)(0, 0).r

    assertEqualsDouble(weight / 4d, 0.5, tolerance, "post-selection probability")

    val normalized = out.scale(1 / Math.sqrt(weight))

    assertEqualsDouble(normalized(0, 0).r, inverseRoot2, tolerance)
    assertEqualsDouble(normalized(1, 0).r, 0d, tolerance)
    assertEqualsDouble(normalized(2, 0).r, 0d, tolerance)
    assertEqualsDouble(normalized(3, 0).r, inverseRoot2, tolerance)
  }

  // -- Allen, the state-discrimination circuit -------------------------------------------------

  /** `U = CH_{CR,CV} SWAP`, which lets a Deutsch curve tell `|0>` from `|->` with certainty. */
  private val discriminator =
    Circuits.inOrder(
      Circuits.swap(qubits = 2, first = 0, second = 1),
      Circuits.controlledHadamard(qubits = 2, control = 0, target = 1)
    )

  test("Allen: a Deutsch curve distinguishes the non-orthogonal states |0> and |-> perfectly") {
    assert(DeutschCtc.solve(discriminator, zero, qubit).output.approximatelyEquals(zero, 1e-6))
    assert(DeutschCtc.solve(discriminator, minus, qubit).output.approximatelyEquals(one, 1e-6))
  }

  test("Allen: the same circuit under post-selection has loop operator |0><0| + |1><-|") {
    val p = PctcLoop.loopOperator(discriminator, Vector(2, 2), ctcSubsystem = 1)
    val expected = ComplexMatrix.outer(
      ComplexMatrix.column(Vector(Complex(1), Complex(0))),
      ComplexMatrix.column(Vector(Complex(1), Complex(0)))
    ).plus(
      ComplexMatrix.outer(
        ComplexMatrix.column(Vector(Complex(0), Complex(1))),
        ComplexMatrix.column(Vector(Complex(inverseRoot2), Complex(-inverseRoot2)))
      )
    )

    assert(p.approximatelyEquals(expected, tolerance), p.toHumanString)
  }

  test("Allen: that loop operator sends |+> to |0> and |1> to |1>, up to phase") {
    val p = PctcLoop.loopOperator(discriminator, Vector(2, 2), ctcSubsystem = 1)

    for (input, expected) <- Seq(plus -> zero, one -> one) do
      val transformed = p.multiply(input.matrix).multiply(p.adjoint)
      val weight = transformed.trace.r

      assert(weight > tolerance)
      assert(DensityMatrix(qubit, transformed.scale(1 / weight)).approximatelyEquals(expected, 1e-9))
  }

  // -- Consistency between the two models ------------------------------------------------------

  test("Deutsch is always satisfiable where post-selection forbids the paradox outright") {
    val u = Circuits.cnot(qubits = 2, control = 0, target = 1)

    // Deutsch: a consistent solution exists.
    assert(DeutschCtc.consistentStates(u, one, qubit).particular.realTrace > 0.5)

    // Post-selected: the loop operator annihilates the time traveller, so the branch never happens.
    val c = PctcLoop.loopOperator(u, Vector(2, 2), ctcSubsystem = 1)
    val transformed = c.multiply(one.matrix).multiply(c.adjoint)

    assertEqualsDouble(transformed.trace.r, 0d, tolerance, "the paradoxical branch has zero weight")
  }

  // -- The side-by-side table ------------------------------------------------------------------

  test("the comparison table reproduces the published contrast on each circuit") {
    val rows = CtcComparison.all(SearchGrid(steps = 13, rounds = 30, radius = 1d))
    val byScenario = rows.map(r => r.scenario -> r).toMap

    assert(rows.nonEmpty)
    assert(rows.forall(r => r.deutschEntropy >= -1e-9), "entropy is never negative")
    assert(rows.forall(r => r.pctcProbability >= -1e-12 && r.pctcProbability <= 1 + 1e-12))

    // Allen: Deutsch leaves the unproven-theorem output classically correlated, purity 1/2, while
    // post-selection leaves it pure and entangled.
    val unproven = byScenario("unproven theorem, |0>_B|0>_M")

    assertEqualsDouble(unproven.deutschOutputPurity, 0.5, 1e-6)
    assertEqualsDouble(unproven.pctcOutputPurity, 1d, 1e-6)
    assertEqualsDouble(unproven.pctcProbability, 0.5, 1e-9)
    assertEquals(unproven.deutschParameters, 1)
    assert(unproven.modelsDisagree)

    // The grandfather paradox: Deutsch always has a consistent history, post-selection has none.
    val paradox = byScenario("grandfather CNOT, traveller |1>")

    assert(paradox.pctcForbids, "post-selection gives the paradoxical branch zero weight")
    assertEqualsDouble(paradox.deutschEntropy, 1d, 1e-6, "Deutsch pays for it with a maximally mixed curve")

    // Ry(pi) is the same paradox written as a rotation.
    assert(byScenario("Ry(theta) on the curve, theta = 1.0 pi").pctcForbids)
    assertEqualsDouble(byScenario("Ry(theta) on the curve, theta = 0.5 pi").pctcProbability, 0.5, 1e-9)
    assertEqualsDouble(byScenario("Ry(theta) on the curve, theta = 0.0 pi").pctcProbability, 1d, 1e-9)
  }

  test("no scenario leaves Deutsch without a consistent history") {
    for row <- CtcComparison.all(SearchGrid(steps = 9, rounds = 20, radius = 1d)) do
      assert(row.deutschEntropy.isFinite, row.scenario)
      assert(row.deutschOutputPurity > 0d, row.scenario)
  }
