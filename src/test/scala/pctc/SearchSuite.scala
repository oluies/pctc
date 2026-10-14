package pctc

import scotty.quantum.density.{ComplexMatrix, Eigen}

class SearchSuite extends munit.FunSuite:
  private val grid = SearchGrid(steps = 9, rounds = 20, radius = 1e3)

  test("a search grid rejects degenerate settings") {
    intercept[IllegalArgumentException](SearchGrid(steps = 1, rounds = 2, radius = 1d))
    intercept[IllegalArgumentException](SearchGrid(steps = 3, rounds = -1, radius = 1d))
    intercept[IllegalArgumentException](SearchGrid(steps = 3, rounds = 2, radius = 0d))
  }

  test("the coordinate search finds the minimum of a smooth two-dimensional objective") {
    val box = Vector((-1d, 1d), (-1d, 1d))
    val objective: CoordinateSearch.Objective = (point, ceiling) =>
      val value = Math.pow(point(0) - 0.3137, 2) + Math.pow(point(1) + 0.7129, 2)

      Option.when(value < ceiling)(value)

    val result = CoordinateSearch.minimize(box, grid, objective)

    assertEqualsDouble(result.value, 0d, 1e-12)
    assertEqualsDouble(result.point.get(0), 0.3137, 1e-6)
    assertEqualsDouble(result.point.get(1), -0.7129, 1e-6)
  }

  test("the coordinate search reports no point when nothing is feasible") {
    val result = CoordinateSearch.minimize(Vector((0d, 1d)), grid, (_, _) => None)

    assertEquals(result.point, None)
    assertEqualsDouble(result.value, OperatorOrder.Ceiling, 0d)
    assert(result.evaluations > 0)
  }

  test("every point of a Bloch state family is a unit-trace state") {
    for dimension <- Seq(2, 3) do
      val family = Bloch.stateFamily(dimension)

      for point <- samples(family) do
        val operator = family.build(point)

        assertEqualsDouble(operator.trace.r, 1d, 1e-12, s"$dimension at $point")
        assert(operator.isHermitian(1e-12), s"$dimension at $point")
        assert(Eigen.minEigenvalue(operator) >= -1e-12, s"$dimension at $point: ${operator.toHumanString}")
  }

  test("every point of a Hermitian family has unit trace and may be indefinite") {
    for dimension <- Seq(2, 3) do
      val family = Bloch.hermitianFamily(dimension, 1e3)
      val points = samples(family)

      for point <- points do
        val operator = family.build(point)

        assertEqualsDouble(operator.trace.r, 1d, 1e-9, s"$dimension at $point")
        assert(operator.isHermitian(1e-9), s"$dimension at $point")

      assert(points.exists(p => Eigen.minEigenvalue(family.build(p)) < -0.5), s"dimension $dimension")
  }

  test("the Hermitian family reaches the Bloch radius it advertises") {
    val family = Bloch.hermitianFamily(2, 1e3)
    val top = family.box.map(_._2)
    val extreme = family.build(Vector(top(0), 0d, 0d))

    // Bloch radius r means eigenvalues (1 +- r) / 2.
    assert(Eigen.maxEigenvalue(extreme) > 400d, extreme.toHumanString)
  }

  test("a block-diagonal qutrit keeps the flag sector separate") {
    val operator = Bloch.blockQutrit(Bloch.qubitOperator(0.25, 0.1, -0.2, 0.05), 0.25)

    assertEqualsDouble(operator(2, 2).r, 0.75, 1e-15)
    assertEqualsDouble(operator(0, 2).r, 0d, 0d)
    assertEqualsDouble(operator(2, 1).r, 0d, 0d)
    assertEqualsDouble(operator.trace.r, 1d, 1e-15)
  }

  private def samples(family: OperatorFamily): Vector[Vector[Double]] =
    val steps = 4

    family.box
      .map((low, high) => Vector.tabulate(steps)(i => low + (high - low) * i / (steps - 1)))
      .foldLeft(Vector(Vector.empty[Double]))((acc, values) => acc.flatMap(p => values.map(p :+ _)))
