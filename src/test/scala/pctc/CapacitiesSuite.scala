package pctc

import munit.FunSuite
import scotty.quantum.density.{Channel, DensityMatrix}

class CapacitiesSuite extends FunSuite:
  private val grid = SearchGrid(steps = 9, rounds = 30, radius = 1e3)
  private val parameters = (1 to 10).map(_ / 10d)

  private def check(name: String, estimated: Double, reference: Double, tolerance: Double): Unit =
    if reference.isInfinite then assert(estimated.isInfinite, s"$name: expected infinity, got $estimated")
    else
      assert(
        Math.abs(estimated - reference) <= tolerance,
        s"$name: estimated $estimated, reference $reference"
      )

  test("the Choi state of the identity channel is the maximally entangled state") {
    assert(Capacities.choiState(Channel.identity(2)).approximatelyEquals(DensityMatrix.maximallyEntangled(2), 1e-12))
  }

  test("I_max matches the closed form for the depolarizing channel") {
    for p <- parameters do
      check(
        s"I_max(D_$p)",
        Capacities.maxInformation(Channel.depolarizing(p), grid).bits,
        ClosedForms.depolarizingMaxInformation(p),
        1e-6
      )
  }

  test("I_doe matches the closed form for the depolarizing channel") {
    for p <- parameters do
      check(
        s"I_doe(D_$p)",
        Capacities.doeblinInformation(Channel.depolarizing(p), grid).bits,
        ClosedForms.depolarizingDoeblinInformation(p),
        1e-6
      )
  }

  test("I_max and I_doe match the closed forms for the amplitude-damping channel") {
    for gamma <- parameters do
      check(
        s"I_max(A_$gamma)",
        Capacities.maxInformation(Channel.amplitudeDamping(gamma), grid).bits,
        ClosedForms.amplitudeDampingMaxInformation(gamma),
        1e-5
      )
      check(
        s"I_doe(A_$gamma)",
        Capacities.doeblinInformation(Channel.amplitudeDamping(gamma), grid).bits,
        ClosedForms.amplitudeDampingDoeblinInformation(gamma),
        1e-5
      )
  }

  test("I_max and I_doe match the closed forms for the erasure channel") {
    for p <- parameters do
      check(
        s"I_max(E_$p)",
        Capacities.maxInformation(Channel.erasure(p), grid).bits,
        ClosedForms.erasureMaxInformation(p),
        1e-5
      )
      check(
        s"I_doe(E_$p)",
        Capacities.doeblinInformation(Channel.erasure(p), grid).bits,
        ClosedForms.erasureDoeblinInformation(p),
        1e-5
      )
  }

  test("a noiseless CTC has infinite Doeblin information and two bits of max-information") {
    check("I_max(id)", Capacities.maxInformation(Channel.identity(2), grid).bits, 2d, 1e-6)
    check("I_doe(id)", Capacities.doeblinInformation(Channel.identity(2), grid).bits, Double.PositiveInfinity, 0d)
  }
