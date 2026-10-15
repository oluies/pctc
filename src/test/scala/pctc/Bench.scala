package pctc

import scotty.quantum.density.{Channel, Cholesky, ComplexMatrix, DensityMatrix, Eigen}

object Bench:
  private def time[A](label: String)(body: => A): A =
    val start = System.nanoTime()
    val result = body
    println(f"$label%-44s ${(System.nanoTime() - start) / 1e9}%9.3f s")
    result

  def main(args: Array[String]): Unit =
    println("--- PSD test: Jacobi eigenvalues vs Cholesky ---")
    for dim <- Seq(4, 6) do
      val m = DensityMatrix.maximallyMixed(Vector(dim)).matrix
      val n = 1000000
      time(f"${dim}x$dim Eigen.minEigenvalue   x $n%,d") {
        var i = 0; var s = 0
        while i < n do { if Eigen.minEigenvalue(m) > -1e-9 then s += 1; i += 1 }; s
      }
      time(f"${dim}x$dim Cholesky              x $n%,d") {
        var i = 0; var s = 0
        while i < n do { if Cholesky.isPositiveDefinite(m, 1e-12) then s += 1; i += 1 }; s
      }

    println("\n--- where an evaluation goes (6x6, erasure) ---")
    val choi = Capacities.choiState(Channel.erasure(0.3)).matrix
    val pi = ComplexMatrix.identity(2).scale(0.5)
    val family = Bloch.hermitianFamily(3, 1e3)
    val point = Vector(0.1, 0.2, 0.3, 0.4)
    val n2 = 1000000
    time(f"family.build                         x $n2%,d") {
      var i = 0; var s = 0.0
      while i < n2 do { s += family.build(point)(0, 0).r; i += 1 }; s
    }
    val operator = family.build(point)
    time(f"pi tensor operator                   x $n2%,d") {
      var i = 0; var s = 0.0
      while i < n2 do { s += pi.tensor(operator)(0, 0).r; i += 1 }; s
    }
    val padded = pi.tensor(operator)
    time(f"OperatorOrder.dominates              x $n2%,d") {
      var i = 0; var s = 0
      while i < n2 do { if OperatorOrder.dominates(padded, choi, 3d) then s += 1; i += 1 }; s
    }
    time(f"Cholesky alone on the 6x6            x $n2%,d") {
      var i = 0; var s = 0
      while i < n2 do { if Cholesky.isPositiveDefinite(choi, 1e-12) then s += 1; i += 1 }; s
    }

    println("\n--- capacity sweep ---")
    for (steps, rounds) <- Seq((9, 30), (13, 40), (17, 50)) do
      val rows = time(s"sweep steps=$steps rounds=$rounds") {
        Demos.retrocausalSweep(Vector.tabulate(11)(_ / 10d), SearchGrid(steps, rounds, 1e3))
      }
      val worst = rows
        .filter(r => r.maxInformationReference.isFinite && r.doeblinInformationReference.isFinite)
        .map(r =>
          Math.max(
            Math.abs(r.maxInformation - r.maxInformationReference),
            Math.abs(r.doeblinInformation - r.doeblinInformationReference)
          )
        )
        .max
      println(f"    evaluations ${rows.map(_.evaluations.toLong).sum}%,d   worst error $worst%.3e")
