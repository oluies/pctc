package pctc

import munit.FunSuite
import scotty.quantum.density.{Channel, ComplexMatrix, DensityMatrix}
import scotty.quantum.gate.StandardGate.{CNOT, H, RY, X}
import scotty.quantum.math.Complex

class PctcLoopSuite extends FunSuite:
  private val tolerance = 1e-10

  private def trivial = DensityMatrix(Vector(1), ComplexMatrix.identity(1))

  test("Eq. (1) with a unitary T and a noiseless CTC matches the Tr_CTC(U) formula") {
    // U acts on (E, B) -> (F, A); subsystem 1 is the chronology-violating one.
    val unitaries = Seq(
      "CNOT" -> ComplexMatrix.fromRows(
        Seq(
          Seq(Complex(1), Complex(0), Complex(0), Complex(0)),
          Seq(Complex(0), Complex(1), Complex(0), Complex(0)),
          Seq(Complex(0), Complex(0), Complex(0), Complex(1)),
          Seq(Complex(0), Complex(0), Complex(1), Complex(0))
        )
      ),
      "H (x) RY" -> ComplexMatrix.fromRaw(H(0).matrix).tensor(ComplexMatrix.fromRaw(RY(0.9, 0).matrix)),
      "swap" -> ComplexMatrix.fromRows(
        Seq(
          Seq(Complex(1), Complex(0), Complex(0), Complex(0)),
          Seq(Complex(0), Complex(0), Complex(1), Complex(0)),
          Seq(Complex(0), Complex(1), Complex(0), Complex(0)),
          Seq(Complex(0), Complex(0), Complex(0), Complex(1))
        )
      )
    )

    val states = Seq(
      DensityMatrix.basis(Vector(2), 0),
      DensityMatrix.basis(Vector(2), 1),
      DensityMatrix.pure(Vector(2), Vector(Complex(0.6), Complex(0, 0.8))),
      DensityMatrix.maximallyMixed(Vector(2))
    )

    for (name, u) <- unitaries; rho <- states do
      val t = PctcLoop.unitaryMap(u, Vector(2, 2), Vector(2, 2))
      val outcome = PctcLoop.evaluate(trivial.tensor(rho), t, Channel.identity(2))

      val c = PctcLoop.loopOperator(u, Vector(2, 2), ctcSubsystem = 1)
      val transformed = c.multiply(rho.matrix).multiply(c.adjoint)
      val weight = transformed.trace.r

      assertEqualsDouble(outcome.postSelectionProbability, weight / 4d, tolerance, name)

      if weight > tolerance then
        assert(
          outcome
            .normalizedOrThrow
            .approximatelyEquals(DensityMatrix(Vector(1, 2), transformed.scale(1 / weight)), tolerance),
          name
        )
      else assertEquals(outcome.normalized(), None, s"$name: a degenerate loop has no normalised output")
  }

  test("the grandfather paradox probability is cos^2(theta/2)") {
    for step <- 0 to 32 do
      val theta = Math.PI * step / 32d
      val t = PctcLoop.unitaryMap(ComplexMatrix.fromRaw(RY(theta, 0).matrix), Vector(1, 2), Vector(1, 2))
      val outcome = PctcLoop.evaluate(trivial.tensor(trivial), t, Channel.identity(2))
      val expected = Math.pow(Math.cos(theta / 2), 2)

      assertEqualsDouble(outcome.postSelectionProbability, expected, 1e-12, s"theta = $theta")
  }

  test("a noiseless loop with T = swap transmits an entangled message perfectly") {
    val swap = ComplexMatrix.fromRows(
      Seq(
        Seq(Complex(1), Complex(0), Complex(0), Complex(0)),
        Seq(Complex(0), Complex(0), Complex(1), Complex(0)),
        Seq(Complex(0), Complex(1), Complex(0), Complex(0)),
        Seq(Complex(0), Complex(0), Complex(0), Complex(1))
      )
    )

    val bell = DensityMatrix.maximallyEntangled(2)
    val t = PctcLoop.unitaryMap(swap, Vector(2, 2), Vector(2, 2))
    val outcome = PctcLoop.evaluate(bell, t, Channel.identity(2))

    assertEqualsDouble(outcome.postSelectionProbability, 0.25, tolerance)
    assert(outcome.normalizedOrThrow.approximatelyEquals(bell, tolerance))
  }
