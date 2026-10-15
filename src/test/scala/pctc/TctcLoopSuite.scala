package pctc

import scotty.quantum.density.{ComplexMatrix, DensityMatrix}
import scotty.quantum.math.Complex

/** Checked against J.-M. A. Allen, Phys. Rev. A 90, 042107 (2014), arXiv:1401.4933: Eq. (33) for the
  * equation of motion and Eq. (35) for the output it gives on the unproven-theorem circuit.
  */
class TctcLoopSuite extends munit.FunSuite:
  private val tolerance = 1e-12
  private val qubit = Vector(2)
  private val root2Inverse = 1 / Math.sqrt(2)

  test("Allen Eq. (35): the unproven-theorem circuit's T-CTC output") {
    val outcome = TctcLoop.evaluate(CtcComparison.unprovenTheorem, DensityMatrix.basis(Vector(2, 2), 0), qubit)
    val got = outcome.normalized

    // rho_f = 1/2 |00><00| + 1/4 |00><11| + 1/4 |11><00| + 1/2 |11><11|
    assertEqualsDouble(got(0, 0).r, 0.5, tolerance)
    assertEqualsDouble(got(3, 3).r, 0.5, tolerance)
    assertEqualsDouble(got(0, 3).r, 0.25, tolerance)
    assertEqualsDouble(got(3, 0).r, 0.25, tolerance)

    for row <- 0 until 4; col <- 0 until 4 if Set(row, col).exists(i => i == 1 || i == 2) do
      assertEqualsDouble(got(row, col).r, 0d, tolerance, s"($row, $col) should vanish")
      assertEqualsDouble(got(row, col).i, 0d, tolerance, s"($row, $col) should vanish")
  }

  test("that output is a genuine state, and strictly between the other two models in purity") {
    val state = DensityMatrix.basis(Vector(2, 2), 0)
    val tctc = TctcLoop.evaluate(CtcComparison.unprovenTheorem, state, qubit).normalized

    assert(tctc.isHermitian(tolerance))
    assert(tctc.isPositiveSemidefinite(1e-12))
    assertEqualsDouble(tctc.realTrace, 1d, tolerance)

    // Allen: Deutsch decoheres it completely, post-selection leaves it pure, T-CTC lands between.
    assertEqualsDouble(CtcComparison.purity(tctc), 0.625, 1e-12)
    assert(CtcComparison.purity(tctc) > 0.5 && CtcComparison.purity(tctc) < 1d)
  }

  test("Eq. (33) splits into the post-selected term and an ordinary channel") {
    val state = DensityMatrix.basis(Vector(2, 2), 0)
    val outcome = TctcLoop.evaluate(CtcComparison.unprovenTheorem, state, qubit)

    // The coherent term is exactly the post-selected numerator, P rho P*.
    val p = PctcLoop.loopOperator(CtcComparison.unprovenTheorem, Vector(2, 2, 2), ctcSubsystem = 2)
    val expected = DensityMatrix(Vector(2, 2), p.multiply(state.matrix).multiply(p.adjoint))

    assert(outcome.coherent.approximatelyEquals(expected, tolerance))
    assert(outcome.incoherent.isPositiveSemidefinite(1e-12), "the second term is a channel output")
    assertEqualsDouble(
      outcome.coherent.realTrace + outcome.incoherent.realTrace,
      outcome.weight,
      tolerance
    )
  }

  test("a T-CTC never forbids a history, where post-selection forbids the paradox") {
    val flip = Circuits.onQubit(qubits = 2, target = 1, Circuits.rotationY(Math.PI))
    val traveller = DensityMatrix.basis(qubit, 0)

    // Post-selected: the branch has zero weight.
    val p = PctcLoop.loopOperator(flip, Vector(2, 2), ctcSubsystem = 1)

    assertEqualsDouble(p.multiply(traveller.matrix).multiply(p.adjoint).trace.r, 0d, tolerance)

    // T-CTC: still a perfectly good state, carried entirely by the channel term.
    val outcome = TctcLoop.evaluate(flip, traveller, qubit)

    assert(outcome.weight > 0d)
    assertEqualsDouble(outcome.coherentFraction, 0d, tolerance, "nothing comes from the P-CTC term")
    assert(outcome.normalized.isPositiveSemidefinite(1e-12))
    assertEqualsDouble(outcome.normalized.realTrace, 1d, tolerance)
  }

  test("a trivial loop reduces to the input, since both terms agree there") {
    for index <- 0 to 1 do
      val state = DensityMatrix.basis(qubit, index)
      val outcome = TctcLoop.evaluate(ComplexMatrix.identity(4), state, qubit)

      assert(outcome.normalized.approximatelyEquals(state, tolerance), s"basis $index")
      assertEqualsDouble(outcome.coherentFraction, 4d / 6d, 1e-12, "P contributes d^2 against d")
  }

  test("T-CTCs cannot distinguish |0> from |-> with certainty, unlike the other two models") {
    val minus = DensityMatrix.pure(qubit, Seq(Complex(root2Inverse), Complex(-root2Inverse)))

    val fromZero = TctcLoop.evaluate(CtcComparison.discriminator, DensityMatrix.basis(qubit, 0), qubit).normalized
    val fromMinus = TctcLoop.evaluate(CtcComparison.discriminator, minus, qubit).normalized

    // Deutsch and post-selection both send these to orthogonal outputs; a T-CTC does not.
    val overlap = fromZero.multiply(fromMinus).realTrace

    assert(overlap > 0.05, s"outputs stay non-orthogonal, Tr[rho_0 rho_-] = $overlap")
    assert(CtcComparison.purity(fromZero) < 1d - 1e-9, "and neither output is pure")
  }
