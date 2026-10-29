package pctc

import org.scalajs.dom
import org.scalajs.dom.{document, html}
import scala.scalajs.js.timers.setTimeout

/** Drives the page: reads the grid controls, runs the demonstrations, writes the tables out.
  *
  * Everything is computed in the browser from the same sources the command-line program uses. The
  * heavy part is the capacity sweep, so the page renders a "computing" state first and starts the
  * work on the next tick, otherwise the browser paints nothing until it finishes.
  */
object WebMain:
  private def element(id: String): html.Element =
    document.getElementById(id).asInstanceOf[html.Element]

  private def input(id: String): html.Input = document.getElementById(id).asInstanceOf[html.Input]

  /** The field's value, held to its `min` and `max`: an out-of-range grid either fails a `require`,
    * which would leave the page saying "computing" for ever, or locks the tab for minutes. The clamped
    * value is written back so the field shows what was actually used.
    */
  private def intValue(id: String, fallback: Int): Int =
    Option(input(id)).fold(fallback) { field =>
      val low = field.min.toIntOption.getOrElse(Int.MinValue)
      val high = field.max.toIntOption.getOrElse(Int.MaxValue)
      val value = field.value.toIntOption.getOrElse(fallback).max(low).min(high)

      field.value = value.toString
      value
    }

  private def write(id: String, markup: String): Unit = element(id).innerHTML = markup

  private def timed[A](body: => A): (A, Double) =
    val start = System.nanoTime()
    val result = body

    (result, (System.nanoTime() - start) / 1e6)

  /** The cheap demonstrations, which are exact and need no grid at all. */
  private def renderExact(): Unit =
    val (grandfather, grandfatherMs) = timed(Demos.grandfather(intValue("theta-steps", 16)))
    val (teleportation, teleportationMs) = timed(Demos.postselectedTeleportation)

    write("grandfather", WebTables.grandfather(grandfather))
    write("grandfather-time", f"$grandfatherMs%.0f ms")
    write("teleportation", WebTables.teleportation(teleportation))
    write("teleportation-time", f"$teleportationMs%.0f ms")

    val (repetition, repetitionMs) = timed(RepetitionCode.rows)

    write("repetition", WebTables.repetition(repetition))
    write("repetition-time", f"$repetitionMs%.0f ms")

  /** The searches, which are the only part that takes real time. */
  private def renderSearches(): Unit =
    val grid = SearchGrid(
      steps = intValue("steps", 7),
      rounds = intValue("rounds", 20),
      radius = 1e3
    )

    write("retrocausal", "<p class=\"working\">computing the capacity searches...</p>")
    write("comparison", "<p class=\"working\">computing the model comparison...</p>")
    write("grid-label", grid.label)
    // The timing lines are live regions, so a screen reader hears the search start and finish.
    write("retrocausal-time", "computing")
    write("comparison-time", "computing")

    setTimeout(0) {
      val parameters = Vector.tabulate(11)(_ / 10d)
      val (rows, sweepMs) = timed(Demos.retrocausalSweep(parameters, grid))

      write("retrocausal", WebTables.retrocausal(rows))
      val points = TextTable.grouped(rows.map(_.evaluations.toLong).sum)

      write("retrocausal-time", f"$sweepMs%.0f ms, $points grid points")

      val worst = rows
        .filter(r => r.maxInformationReference.isFinite && r.doeblinInformationReference.isFinite)
        .map(r =>
          Math.max(
            Math.abs(r.maxInformation - r.maxInformationReference),
            Math.abs(r.doeblinInformation - r.doeblinInformationReference)
          )
        )
        .max

      write("retrocausal-error", TextTable.scientific(worst))

      setTimeout(0) {
        val (comparison, comparisonMs) = timed(CtcComparison.all(SearchGrid(grid.steps, grid.rounds, 1d)))

        write("comparison", WebTables.comparison(comparison))
        write("comparison-time", f"$comparisonMs%.0f ms")
      }
    }

  def main(args: Array[String]): Unit =
    renderExact()
    renderSearches()

    Option(document.getElementById("run")).foreach { button =>
      button.addEventListener("click", (_: dom.Event) => { renderExact(); renderSearches() })
    }
