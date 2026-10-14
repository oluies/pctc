package pctc

class DemosSuite extends munit.FunSuite:
  private val grid = SearchGrid(steps = 7, rounds = 20, radius = 1e3)

  test("the grandfather sweep matches cos^2(theta/2) to 1e-9") {
    val rows = Demos.grandfather(60)

    assertEquals(rows.length, 61)
    assertEqualsDouble(rows.head.theta, 0d, 0d)
    assertEqualsDouble(rows.last.theta, Math.PI, 1e-15)

    for row <- rows do
      assertEqualsDouble(
        row.postSelectionProbability,
        Math.pow(Math.cos(row.theta / 2), 2),
        1e-9,
        s"theta = ${row.theta}"
      )
      assert(row.absoluteError <= 1e-9, s"theta = ${row.theta}, error ${row.absoluteError}")
  }

  test("the grandfather sweep starts at one and vanishes at theta = pi") {
    val rows = Demos.grandfather(20)

    assertEqualsDouble(rows.head.postSelectionProbability, 1d, 1e-12)
    assertEqualsDouble(rows.last.postSelectionProbability, 0d, 1e-12)
  }

  test("the post-selection probability decreases monotonically in theta") {
    val probabilities = Demos.grandfather(40).map(_.postSelectionProbability)

    assert(probabilities.lazyZip(probabilities.tail).forall(_ > _ - 1e-15))
  }

  test("teleportation succeeds with probability one quarter before renormalisation") {
    for row <- Demos.postselectedTeleportation do
      assertEqualsDouble(row.postSelectionProbability, 0.25, 1e-12, row.stage)
  }

  test("the loop renormalises the teleportation probability to one and the message survives") {
    val rows = Demos.postselectedTeleportation

    assertEquals(rows.length, 4)

    for row <- rows do
      assertEqualsDouble(row.messageFidelity, 1d, 1e-12, row.stage)

      if row.stage.startsWith("P-CTC loop") then
        assertEqualsDouble(row.renormalizedProbability, 1d, 1e-15, row.stage)
      else assertEqualsDouble(row.renormalizedProbability, 0.25, 1e-12, row.stage)
  }

  test("every retrocausal row agrees with the paper's closed forms") {
    val rows = Demos.retrocausalSweep(Vector.tabulate(11)(_ / 10d), grid)

    assertEquals(rows.length, 33)

    for row <- rows do
      val label = s"${row.channel} ${row.parameterName} = ${row.parameter}"

      if row.maxInformationReference.isInfinite then assert(row.maxInformation.isInfinite, label)
      else assertEqualsDouble(row.maxInformation, row.maxInformationReference, 1e-4, s"I_max, $label")

      if row.doeblinInformationReference.isInfinite then assert(row.doeblinInformation.isInfinite, label)
      else assertEqualsDouble(row.doeblinInformation, row.doeblinInformationReference, 1e-3, s"I_doe, $label")
  }

  test("a noiseless CTC gives an infinite retrocausal capacity and a quarter post-selection") {
    val rows = Demos.retrocausalSweep(Vector(0d), grid)

    for row <- rows do
      assert(row.retrocausalCapacity.isInfinite, row.channel)
      assertEqualsDouble(row.postSelectionProbability, 0.25, 1e-6, row.channel)
      assert(row.exceedsOrdinaryCapacity, row.channel)
  }

  test("a fully noisy CTC carries nothing and ties the ordinary capacity") {
    for spec <- ChannelSpec.all do
      val row = Demos.retrocausal(spec, 1d, grid)

      assertEqualsDouble(row.retrocausalCapacity, 0d, 1e-6, spec.name)
      assertEqualsDouble(row.ordinaryClassicalCapacity, 0d, 1e-6, spec.name)
      assert(!row.exceedsOrdinaryCapacity, spec.name)
  }
