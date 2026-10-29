package pctc

class RepetitionCodeSuite extends munit.FunSuite:
  private def assertRelative(actual: Double, expected: Double, relative: Double)(using munit.Location): Unit =
    assert(
      Math.abs(actual - expected) <= relative * Math.abs(expected) + 1e-20,
      s"$actual differs from $expected by more than $relative relative"
    )

  // 1e-5 relative because of the smallest row: at n = 7 and p = 5e-4 the failure rate is 2.2e-12, a
  // sum of diagonal entries near 1e-12 in a matrix whose entries are near 1e-2, and double rounding
  // leaves it accurate to about 7.5e-6. Every other row agrees to better than 1e-7.
  test("every simulated row matches its closed form") {
    RepetitionCode.rows.foreach: row =>
      assertRelative(row.simulated, row.closedForm, 1e-5)
  }

  test("a single qubit fails with exactly the flip probability under either noise") {
    for p <- Seq(0d, 0.1, 0.5, 1d) do
      assertEqualsDouble(RepetitionCode.phaseFlipFailure(1, p), p, 1e-12)
      assertEqualsDouble(RepetitionCode.bitFlipFailure(1, p), p, 1e-12)
  }

  test("distance suppresses phase flips below threshold and amplifies them above it") {
    def failures(p: Double) = Seq(1, 3, 5, 7).map(RepetitionCode.phaseFlipFailure(_, p))

    assert(failures(0.1).sliding(2).forall { case Seq(a, b) => b < a })
    assert(failures(0.6).sliding(2).forall { case Seq(a, b) => b > a })
    failures(0.5).foreach(assertEqualsDouble(_, 0.5, 1e-12))
  }

  test("at the paper's p_Z, each step of distance buys about three orders of magnitude") {
    val d5 = RepetitionCode.phaseFlipFailure(5, RepetitionCode.PaperPhaseFlip)
    val d7 = RepetitionCode.phaseFlipFailure(7, RepetitionCode.PaperPhaseFlip)

    // 10 p^3 and 35 p^4 to leading order.
    assertRelative(d5, 10 * Math.pow(5e-4, 3), 1e-2)
    assertRelative(d7, 35 * Math.pow(5e-4, 4), 1e-2)
  }

  test("bit flips get worse with distance: about d p, because one flip is a logical error") {
    val p = RepetitionCode.PaperPhaseFlip
    val failures = Seq(1, 3, 5, 7).map(RepetitionCode.bitFlipFailure(_, p))

    assert(failures.sliding(2).forall { case Seq(a, b) => b > a })
    Seq(1, 3, 5, 7).zip(failures).foreach((d, f) => assertRelative(f, d * p, 1e-2))
  }

  test("the code rejects an even distance") {
    intercept[IllegalArgumentException](RepetitionCode.phaseFlipFailure(4, 0.1))
    intercept[IllegalArgumentException](RepetitionCode.bitFlipFailureClosedForm(2, 0.1))
  }
