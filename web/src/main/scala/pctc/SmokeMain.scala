package pctc

import scotty.quantum.density.{Channel, DensityMatrix}

/** Verifies the numerics in a JavaScript runtime.
  *
  * Scala.js cannot read classpath resources and has no DuckDB, so `ReferenceSuite` and
  * `ResultsSuite` cannot run there and the rest of the suite is not worth half-porting. Instead this
  * checks the values the JVM suite pins hardest -- the closed forms of arXiv:2509.08965 and the
  * published outputs of arXiv:1401.4933 -- directly against the linked JavaScript, and exits
  * non-zero on any mismatch so continuous integration can gate on it.
  *
  * The arithmetic is all `Double`, which is IEEE-754 in both runtimes, but `Math.sqrt` and friends
  * map to the host's library rather than Java's, so the last bit or two can differ. The tolerances
  * below are loose enough to allow that and tight enough to catch a real divergence.
  */
object SmokeMain:
  private var failures = 0

  private def check(name: String, actual: Double, expected: Double, tolerance: Double): Unit =
    val ok = (actual.isInfinite && expected.isInfinite) || Math.abs(actual - expected) <= tolerance

    if !ok then failures += 1

    println(f"  ${if ok then "ok  " else "FAIL"}  $name%-46s $actual%18.12f  expected $expected%18.12f")

  def main(args: Array[String]): Unit =
    println("P-CTC numerics, running on JavaScript\n")

    for step <- 0 to 8 do
      val theta = Math.PI * step / 8
      val row = Demos.grandfather(8)(step)

      check(f"grandfather, theta = ${step / 8d}%.3f pi", row.postSelectionProbability, Math.pow(Math.cos(theta / 2), 2), 1e-12)

    for row <- Demos.postselectedTeleportation do
      check(s"teleportation, ${row.stage.take(30)}", row.postSelectionProbability, 0.25, 1e-12)

    val grid = SearchGrid(steps = 9, rounds = 30, radius = 1e3)

    for (name, channel, reference) <- Seq(
        ("I_max depolarizing(0.5)", Channel.depolarizing(0.5), ClosedForms.depolarizingMaxInformation(0.5)),
        ("I_max amplitudeDamping(0.5)", Channel.amplitudeDamping(0.5), ClosedForms.amplitudeDampingMaxInformation(0.5)),
        ("I_max erasure(0.5)", Channel.erasure(0.5), ClosedForms.erasureMaxInformation(0.5))
      )
    do check(name, Capacities.maxInformation(channel, grid).bits, reference, 1e-6)

    for (name, channel, reference) <- Seq(
        ("I_doe depolarizing(0.5)", Channel.depolarizing(0.5), ClosedForms.depolarizingDoeblinInformation(0.5)),
        ("I_doe amplitudeDamping(0.5)", Channel.amplitudeDamping(0.5), ClosedForms.amplitudeDampingDoeblinInformation(0.5)),
        ("I_doe erasure(0.5)", Channel.erasure(0.5), ClosedForms.erasureDoeblinInformation(0.5))
      )
    do check(name, Capacities.doeblinInformation(channel, grid).bits, reference, 1e-5)

    check("I_doe identity is infinite", Capacities.doeblinInformation(Channel.identity(2), grid).bits, Double.PositiveInfinity, 0d)

    // Allen arXiv:1401.4933 Eq. (35) and the Deutsch maximum-entropy solution for the same circuit.
    val unproven = DensityMatrix.basis(Vector(2, 2), 0)

    check(
      "Allen Eq. (35) T-CTC purity",
      CtcComparison.purity(TctcLoop.evaluate(CtcComparison.unprovenTheorem, unproven, Vector(2)).normalized),
      0.625,
      1e-12
    )
    check(
      "Deutsch output purity, same circuit",
      CtcComparison.purity(DeutschCtc.solve(CtcComparison.unprovenTheorem, unproven, Vector(2)).output),
      0.5,
      1e-6
    )

    // The repetition-code side note: the density matrix against its closed forms, 1e-5 relative,
    // as in RepetitionCodeSuite, with an absolute floor that only matters for an exact zero.
    for row <- RepetitionCode.rows do
      check(s"${row.noise}, n = ${row.distance}, p = ${row.probability}", row.simulated, row.closedForm, 1e-5 * row.closedForm + 1e-20)

    println(if failures == 0 then "\nall checks passed" else s"\n$failures check(s) FAILED")

    if failures > 0 then throw new RuntimeException(s"$failures smoke check(s) failed")
