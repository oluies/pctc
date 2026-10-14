package pctc

import scotty.quantum.density.{ComplexMatrix, DensityMatrix}

class OperatorOrderSuite extends munit.FunSuite:
  private val identity = ComplexMatrix.identity(4)
  private val bell = DensityMatrix.maximallyEntangled(2).matrix
  private val mixed = DensityMatrix.maximallyMixed(Vector(2, 2)).matrix

  test("an operator dominates itself at lambda = 1") {
    assert(OperatorOrder.dominates(bell, bell, 0d))
    assertEqualsDouble(OperatorOrder.minLogLambda(bell, bell).get, 0d, 1e-9)
  }

  test("the maximally mixed state needs lambda = 4 to dominate a Bell state") {
    // Phi has largest eigenvalue 1 and I/4 has 1/4, so lambda = 4 and log2 lambda = 2.
    assertEqualsDouble(OperatorOrder.minLogLambda(bell, mixed).get, 2d, 1e-8)
  }

  test("a Bell state cannot dominate the maximally mixed state at any lambda") {
    assertEquals(OperatorOrder.minLogLambda(mixed, bell), None)
  }

  test("scaling the reference shifts log2 lambda by the same number of bits") {
    val base = OperatorOrder.minLogLambda(bell, mixed).get

    assertEqualsDouble(OperatorOrder.minLogLambda(bell, mixed.scale(8d)).get, base - 3d, 1e-8)
    assertEqualsDouble(OperatorOrder.minLogLambda(bell, mixed.scale(0.125)).get, base + 3d, 1e-8)
  }

  test("a negative exponent is handled, not clamped at zero") {
    assertEqualsDouble(OperatorOrder.minLogLambda(mixed.scale(0.25), mixed).get, -2d, 1e-8)
  }

  test("the ceiling doubles as an incumbent and rejects anything that cannot beat it") {
    assertEquals(OperatorOrder.minLogLambda(bell, mixed, ceiling = 1d), None)
    assert(OperatorOrder.minLogLambda(bell, mixed, ceiling = 3d).isDefined)
  }

  test("the positivity tolerance stays well below the ceiling's resolution") {
    assert(OperatorOrder.PositivityTolerance < Math.pow(2d, -OperatorOrder.Ceiling) / 16)
  }

  test("an identity reference needs the largest eigenvalue of the candidate") {
    val candidate = DensityMatrix.maximallyEntangled(2).scale(8d).matrix

    assertEqualsDouble(OperatorOrder.minLogLambda(candidate, identity).get, 3d, 1e-8)
  }
