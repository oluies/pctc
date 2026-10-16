package pctc

/** How much of a run is JIT warm-up? That decides whether ahead-of-time compilation could help. */
object WarmupBench:
  def main(args: Array[String]): Unit =
    val parameters = Vector.tabulate(11)(_ / 10d)
    val grid = SearchGrid.default

    println("  run    time     vs steady state")
    val times = (1 to 8).map { _ =>
      val start = System.nanoTime()
      Demos.retrocausalSweep(parameters, grid)
      (System.nanoTime() - start) / 1e9
    }
    val steady = times.drop(4).sum / 4
    times.zipWithIndex.foreach { (t, i) =>
      println(f"  ${i + 1}%3d  ${t}%6.3f s   ${t / steady}%5.2fx")
    }
    println(f"\n  steady state ${steady}%.3f s; first run carries ${(times.head - steady)}%.3f s of warm-up")
    println(f"  a single sweep is ${(times.head - steady) / times.head * 100}%.0f%% warm-up")
