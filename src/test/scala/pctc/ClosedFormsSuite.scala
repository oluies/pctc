package pctc

class ClosedFormsSuite extends munit.FunSuite:
  test("the binary entropy is zero at the ends and one in the middle") {
    assertEqualsDouble(ClosedForms.binaryEntropy(0d), 0d, 0d)
    assertEqualsDouble(ClosedForms.binaryEntropy(1d), 0d, 0d)
    assertEqualsDouble(ClosedForms.binaryEntropy(0.5), 1d, 1e-15)
    assertEqualsDouble(ClosedForms.binaryEntropy(0.11), ClosedForms.binaryEntropy(0.89), 1e-15)
  }

  test("the max-informations span two bits down to zero") {
    assertEqualsDouble(ClosedForms.depolarizingMaxInformation(0d), 2d, 1e-15)
    assertEqualsDouble(ClosedForms.depolarizingMaxInformation(1d), 0d, 1e-15)
    assertEqualsDouble(ClosedForms.erasureMaxInformation(0d), 2d, 1e-15)
    assertEqualsDouble(ClosedForms.erasureMaxInformation(1d), 0d, 1e-15)
    assertEqualsDouble(ClosedForms.amplitudeDampingMaxInformation(0d), 2d, 1e-15)
    assertEqualsDouble(ClosedForms.amplitudeDampingMaxInformation(1d), 0d, 1e-15)
  }

  test("the Doeblin informations diverge for a noiseless channel and vanish for a useless one") {
    assert(ClosedForms.depolarizingDoeblinInformation(0d).isPosInfinity)
    assert(ClosedForms.erasureDoeblinInformation(0d).isPosInfinity)
    assert(ClosedForms.amplitudeDampingDoeblinInformation(0d).isPosInfinity)

    assertEqualsDouble(ClosedForms.depolarizingDoeblinInformation(1d), 0d, 1e-15)
    assertEqualsDouble(ClosedForms.erasureDoeblinInformation(1d), 0d, 1e-15)
    assertEqualsDouble(ClosedForms.amplitudeDampingDoeblinInformation(1d), 0d, 1e-15)
  }

  test("the ordinary classical capacities run from one bit to zero") {
    assertEqualsDouble(ClosedForms.depolarizingClassicalCapacity(0d), 1d, 1e-15)
    assertEqualsDouble(ClosedForms.depolarizingClassicalCapacity(1d), 0d, 1e-15)
    assertEqualsDouble(ClosedForms.erasureClassicalCapacity(0d), 1d, 1e-15)
    assertEqualsDouble(ClosedForms.erasureClassicalCapacity(1d), 0d, 1e-15)
    assertEqualsDouble(ClosedForms.amplitudeDampingHolevoInformation(0d), 1d, 1e-9)
    assertEqualsDouble(ClosedForms.amplitudeDampingHolevoInformation(1d), 0d, 1e-9)
  }

  test("every capacity is non-increasing in the noise parameter") {
    val parameters = Vector.tabulate(21)(_ / 20d)

    def nonIncreasing(name: String, f: Double => Double): Unit =
      val values = parameters.map(f)

      assert(values.lazyZip(values.tail).forall(_ >= _ - 1e-9), s"$name: $values")

    nonIncreasing("depolarizing C", ClosedForms.depolarizingClassicalCapacity)
    nonIncreasing("erasure C", ClosedForms.erasureClassicalCapacity)
    nonIncreasing("damping chi", ClosedForms.amplitudeDampingHolevoInformation(_))
    nonIncreasing("depolarizing I_max", ClosedForms.depolarizingMaxInformation)
    nonIncreasing("erasure I_max", ClosedForms.erasureMaxInformation)
    nonIncreasing("damping I_max", ClosedForms.amplitudeDampingMaxInformation)
  }

  test("the depolarizing capacity has its textbook value at p = 1/2") {
    // 1 - h_2(1/4) = 1 - 0.811278...
    assertEqualsDouble(ClosedForms.depolarizingClassicalCapacity(0.5), 1 - ClosedForms.binaryEntropy(0.25), 1e-15)
    assertEqualsDouble(ClosedForms.depolarizingClassicalCapacity(0.5), 0.18872187554086717, 1e-12)
  }
