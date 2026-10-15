package pctc

import scala.io.Source
import scotty.quantum.density.Channel

/** Cross-validation against an independent implementation.
  *
  * The rest of the suite checks the numerics against closed forms transcribed from the same paper
  * this program implements, so a transcription slip would pass on both sides. The fixtures read here
  * come from two independent implementations that build the three channels from their Kraus
  * definitions and solve Eqs. (9) and (10) as semidefinite programs:
  *
  *   - `reference/sdp_reference.py`, numpy and cvxpy, solved by CLARABEL (interior point);
  *   - `reference/sdp_reference.jl`, Convex.jl, solved by SCS (first order).
  *
  * Neither shares code, language, modelling layer or solver with the Scala search, or with each
  * other.
  *
  * It is also unrestricted: the solver optimises over the whole Hermitian cone, where the Scala
  * search sweeps a Bloch parametrisation and, for the erasure channel, only block-diagonal
  * operators. Agreement is therefore evidence that neither restriction loses anything.
  *
  * Regenerate with `uv run reference/sdp_reference.py` and
  * `julia --project=reference reference/sdp_reference.jl`. The fixtures are committed so the suite
  * needs neither Python nor Julia.
  */
class ReferenceSuite extends munit.FunSuite:
  private val grid = SearchGrid(steps = 9, rounds = 30, radius = 1e3)

  /** The solvers' own accuracy on these programs, both looser than the Scala search's. */
  private val solverTolerance = 1e-4

  /** How closely the two independent solvers have to agree with each other. */
  private val crossSolverTolerance = 1e-3

  private final case class Solver(name: String, resource: String)

  private val solvers = Vector(
    Solver("cvxpy / CLARABEL", "reference-capacities.tsv"),
    Solver("Convex.jl / SCS", "reference-capacities-julia.tsv")
  )

  /** The numpy reference sweeps the Deutsch family on a fixed 81-point grid where the Scala side
    * refines, so its entropies and purities carry about 1e-4 of grid error. The Scala values are the
    * tighter ones; this bound is the reference's slack, not ours.
    */
  private val modelTolerance = 1e-3

  private final case class Model(
      scenario: String,
      parameters: Int,
      entropy: Double,
      deutschPurity: Double,
      pctcProbability: Double,
      pctcPurity: Double
  )

  private lazy val models: Vector[Model] =
    read("reference-ctc-models.tsv").map(c =>
      Model(c(0), c(1).toInt, number(c(2)), number(c(3)), number(c(4)), if c(5) == "nan" then Double.NaN else number(c(5)))
    )

  private final case class Capacity(channel: String, parameter: Double, iMax: Double, iDoe: Double, defect: Double)

  private final case class ChoiEntry(channel: String, parameter: Double, row: Int, col: Int, re: Double, im: Double)

  private def read(resource: String): Vector[Vector[String]] =
    val source = Source.fromResource(resource)

    try source.getLines().drop(1).filter(_.nonEmpty).map(_.split('\t').toVector).toVector
    finally source.close()

  private def number(text: String): Double = text match
    case "inf"  => Double.PositiveInfinity
    case "-inf" => Double.NegativeInfinity
    case other  => other.toDouble

  private def capacitiesOf(resource: String): Vector[Capacity] =
    read(resource).map(c => Capacity(c(0), number(c(1)), number(c(2)), number(c(3)), number(c(4))))

  private lazy val capacities: Vector[Capacity] = capacitiesOf(solvers.head.resource)

  private lazy val chois: Vector[ChoiEntry] =
    read("reference-choi.tsv").map(c => ChoiEntry(c(0), number(c(1)), c(2).toInt, c(3).toInt, number(c(4)), number(c(5))))

  private def specFor(name: String): ChannelSpec =
    ChannelSpec.all.find(_.name == name).getOrElse(fail(s"no channel named $name"))

  test("the fixtures cover every channel at every parameter") {
    for solver <- solvers do
      val rows = capacitiesOf(solver.resource)

      assertEquals(rows.length, 33, solver.name)
      assertEquals(rows.map(_.channel).distinct.sorted, ChannelSpec.all.map(_.name).sorted, solver.name)
      assertEquals(rows.map(_.parameter).distinct.sorted, Vector.tabulate(11)(_ / 10d), solver.name)

    assert(chois.nonEmpty)
  }

  test("the reference channels satisfy Kraus completeness, as ours do") {
    for solver <- solvers; row <- capacitiesOf(solver.resource) do
      assert(row.defect < 1e-12, s"${solver.name}, ${row.channel} ${row.parameter}: ${row.defect}")
  }

  test("the two independent solvers agree with each other") {
    val Vector(first, second) = solvers.map(s => capacitiesOf(s.resource))

    for (a, b) <- first.zip(second) do
      assertEquals(a.channel -> a.parameter, b.channel -> b.parameter)

      val label = s"${a.channel} = ${a.parameter}"

      assertEquals(a.iMax.isInfinite, b.iMax.isInfinite, s"I_max finiteness, $label")
      assertEquals(a.iDoe.isInfinite, b.iDoe.isInfinite, s"I_doe finiteness, $label")

      if a.iMax.isFinite then assertEqualsDouble(a.iMax, b.iMax, crossSolverTolerance, s"I_max, $label")
      if a.iDoe.isFinite then assertEqualsDouble(a.iDoe, b.iDoe, crossSolverTolerance, s"I_doe, $label")
  }

  test("our Choi operators match numpy's entry by entry") {
    val built = capacities.map(row => (row.channel, row.parameter) -> Capacities.choiState(specFor(row.channel).build(row.parameter))).toMap

    for entry <- chois do
      val actual = built((entry.channel, entry.parameter))(entry.row, entry.col)
      val label = s"${entry.channel} ${entry.parameter} at (${entry.row}, ${entry.col})"

      assertEqualsDouble(actual.r, entry.re, 1e-12, s"$label, real part")
      assertEqualsDouble(actual.i, entry.im, 1e-12, s"$label, imaginary part")
  }

  test("our Choi operators are zero wherever numpy's fixture omits an entry") {
    val present = chois.map(e => (e.channel, e.parameter, e.row, e.col)).toSet

    for row <- capacities do
      val choi = Capacities.choiState(specFor(row.channel).build(row.parameter))

      for r <- 0 until choi.dimension; c <- 0 until choi.dimension
      if !present.contains((row.channel, row.parameter, r, c))
      do
        val entry = choi(r, c)

        assert(
          Math.abs(entry.r) < 1e-15 && Math.abs(entry.i) < 1e-15,
          s"${row.channel} ${row.parameter} at ($r, $c) should vanish but is $entry"
        )
  }

  test("I_max agrees with both independent semidefinite solutions") {
    for solver <- solvers; row <- capacitiesOf(solver.resource) do
      val ours = Capacities.maxInformation(specFor(row.channel).build(row.parameter), grid).bits
      val label = s"I_max, ${solver.name}, ${row.channel} = ${row.parameter}"

      if row.iMax.isInfinite then assert(ours.isInfinite, label)
      else assertEqualsDouble(ours, row.iMax, solverTolerance, label)
  }

  test("I_doe agrees with both independent semidefinite solutions, infinities included") {
    for solver <- solvers; row <- capacitiesOf(solver.resource) do
      val ours = Capacities.doeblinInformation(specFor(row.channel).build(row.parameter), grid).bits
      val label = s"I_doe, ${solver.name}, ${row.channel} = ${row.parameter}"

      if row.iDoe.isInfinite then assert(ours.isInfinite, s"$label: expected infinity, got $ours")
      else assertEqualsDouble(ours, row.iDoe, solverTolerance, label)
  }

  /** Slack for the comparison below.
    *
    * Not zero: a solver sometimes lands exactly on the closed form, where our search carries the
    * bisection tolerance of about 1e-9. The claim worth making is that our search is never the
    * *loose* one by a margin that matters, which this bound catches while a regression to, say,
    * 1e-6 would fail it.
    */
  private val accuracySlack = 1e-8

  test("our search tracks the closed forms at least as tightly as either solver") {
    for solver <- solvers; row <- capacitiesOf(solver.resource) do
      val spec = specFor(row.channel)
      val referenceMax = spec.referenceMaxInformation(row.parameter)
      val referenceDoe = spec.referenceDoeblinInformation(row.parameter)

      if referenceMax.isFinite then
        val ours = Math.abs(Capacities.maxInformation(spec.build(row.parameter), grid).bits - referenceMax)

        assert(
          ours <= Math.abs(row.iMax - referenceMax) + accuracySlack,
          s"I_max, ${solver.name}, ${row.channel} = ${row.parameter}: ours $ours"
        )

      if referenceDoe.isFinite then
        val ours = Math.abs(Capacities.doeblinInformation(spec.build(row.parameter), grid).bits - referenceDoe)

        assert(
          ours <= Math.abs(row.iDoe - referenceDoe) + accuracySlack,
          s"I_doe, ${solver.name}, ${row.channel} = ${row.parameter}: ours $ours"
        )
  }

  // -- Both closed-timelike-curve models, against an independent numpy implementation -----------

  test("the Deutsch and post-selected models agree with numpy on every published circuit") {
    val computed = CtcComparison.all(SearchGrid(steps = 13, rounds = 30, radius = 1d))
      .map(r => r.scenario -> r)
      .toMap

    assertEquals(models.length, 10)
    assertEquals(models.map(_.scenario).toSet, computed.keySet)

    for reference <- models do
      val ours = computed(reference.scenario)
      val label = reference.scenario

      assertEquals(ours.deutschParameters, reference.parameters, s"solution dimension, $label")
      assertEqualsDouble(ours.deutschEntropy, reference.entropy, modelTolerance, s"S(tau), $label")
      assertEqualsDouble(ours.deutschOutputPurity, reference.deutschPurity, modelTolerance, s"Deutsch purity, $label")
      assertEqualsDouble(ours.pctcProbability, reference.pctcProbability, 1e-9, s"P(post-select), $label")

      if reference.pctcPurity.isNaN then
        assert(ours.pctcOutputPurity.isNaN, s"$label: post-selection should forbid this branch")
        assert(ours.pctcForbids, label)
      else assertEqualsDouble(ours.pctcOutputPurity, reference.pctcPurity, modelTolerance, s"P-CTC purity, $label")
  }

  test("the two models genuinely disagree somewhere, so the comparison is not vacuous") {
    val rows = CtcComparison.all(SearchGrid(steps = 9, rounds = 20, radius = 1d))

    assert(rows.exists(_.pctcForbids), "post-selection forbids at least one history")
    assert(rows.forall(_.deutschEntropy.isFinite), "Deutsch always admits one")
    assert(rows.exists(r => Math.abs(r.deutschOutputPurity - r.pctcOutputPurity) > 0.4),
      "and they differ sharply on at least one circuit")
  }
