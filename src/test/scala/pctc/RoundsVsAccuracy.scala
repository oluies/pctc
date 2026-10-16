package pctc

object RoundsVsAccuracy:
  def main(args: Array[String]): Unit =
    val parameters = Vector.tabulate(11)(_ / 10d)

    def worstError(rows: Vector[RetrocausalRow]): Double =
      rows
        .filter(r => r.maxInformationReference.isFinite && r.doeblinInformationReference.isFinite)
        .map(r =>
          Math.max(
            Math.abs(r.maxInformation - r.maxInformationReference),
            Math.abs(r.doeblinInformation - r.doeblinInformationReference)
          )
        )
        .max

    // Warm up the JIT so the first row is not penalised.
    Demos.retrocausalSweep(parameters, SearchGrid(17, 5, 1e3))

    println("  steps  rounds      time   worst error vs closed forms   evaluations")
    for (steps, rounds) <- Seq((9, 20), (9, 30), (9, 40), (11, 25), (11, 30), (13, 25), (13, 30), (17, 25), (17, 30)) do
      val grid = SearchGrid(steps, rounds, 1e3)
      val start = System.nanoTime()
      val rows = Demos.retrocausalSweep(parameters, grid)
      val wall = (System.nanoTime() - start) / 1e9

      println(f"  $steps%5d  $rounds%6d  ${wall}%6.3f s   ${worstError(rows)}%.3e"
        + f"                  ${rows.map(_.evaluations.toLong).sum}%,12d")
