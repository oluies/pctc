package pctc

import scotty.quantum.density.{Cholesky, ComplexMatrix}
import scotty.quantum.math.Complex

/** How much does the boxed `Vector[Complex]` representation cost the inner loop? */
object ReprBench:
  /** The same Cholesky, reading two interleaved primitive arrays instead of ComplexMatrix. */
  private def primitiveFactorizes(
      left: Array[Double],
      right: Array[Double],
      n: Int,
      factor: Double,
      shift: Double
  ): Boolean =
    val re = new Array[Double](n * n)
    val im = new Array[Double](n * n)
    var row = 0
    var definite = true

    while row < n && definite do
      val d = 2 * (row * n + row)
      var pivotSquared = left(d) - factor * right(d) + shift
      var k = 0
      while k < row do
        val lr = re(row * n + k); val li = im(row * n + k)
        pivotSquared -= lr * lr + li * li
        k += 1

      if pivotSquared <= 0d then definite = false
      else
        val pivot = Math.sqrt(pivotSquared)
        re(row * n + row) = pivot
        var below = row + 1
        while below < n do
          val b = 2 * (below * n + row)
          var sumR = left(b) - factor * right(b)
          var sumI = left(b + 1) - factor * right(b + 1)
          var j = 0
          while j < row do
            val br = re(below * n + j); val bi = im(below * n + j)
            val rr = re(row * n + j); val ri = im(row * n + j)
            sumR -= br * rr + bi * ri
            sumI -= bi * rr - br * ri
            j += 1
          re(below * n + row) = sumR / pivot
          im(below * n + row) = sumI / pivot
          below += 1
        row += 1

    definite

  private def interleaved(m: ComplexMatrix): Array[Double] =
    Array.tabulate(m.rows * m.cols * 2) { k =>
      val c = m(k / 2 / m.cols, (k / 2) % m.cols)
      if k % 2 == 0 then c.r else c.i
    }

  def main(args: Array[String]): Unit =
    val choi = Capacities.choiState(scotty.quantum.density.Channel.erasure(0.3)).matrix
    val pi = ComplexMatrix.identity(2).scale(0.5)
    val padded = pi.tensor(Bloch.hermitianFamily(3, 1e3).build(Vector(0.1, 0.2, 0.3, 0.4)))
    val choiFlat = interleaved(choi)
    val paddedFlat = interleaved(padded)
    val n = 2000000
    val factor = Math.pow(2d, -3d)

    def time(label: String)(body: => Boolean): Double =
      var i = 0; var s = 0
      val start = System.nanoTime()
      while i < n do { if body then s += 1; i += 1 }
      val per = (System.nanoTime() - start).toDouble / n / 1000
      println(f"  $label%-50s ${per}%7.3f us/call  (sink $s)")
      per

    println(s"6x6 positive-definiteness, $n calls each:\n")
    val boxed = time("ComplexMatrix-backed (what ships)") {
      Cholesky.differenceIsPositiveDefinite(choi, padded, factor, 1e-12)
    }
    val primitive = time("primitive Array[Double], same algorithm") {
      primitiveFactorizes(choiFlat, paddedFlat, 6, factor, 1e-12 * 0.5)
    }
    println(f"\n  primitive is ${boxed / primitive}%.2fx faster than the boxed representation")
