package pctc

import scala.concurrent.duration.Duration
import scala.concurrent.{Await, ExecutionContext, Future}

/** The one thing this program does differently per platform.
  *
  * Splitting the coarse grid across threads is worth about 1.2x on a whole sweep, measured back to
  * back with a control, so it is worth keeping on the JVM. A browser has no `Await`, so the Scala.js
  * build supplies its own sequential version of this object; see `web/src/main/scala`.
  */
object Platform:
  val threads: Int = Runtime.getRuntime.availableProcessors

  /** Runs `scan` over each range, in parallel where the platform allows it. */
  def mapRanges[A](ranges: Vector[(Int, Int)])(scan: (Int, Int) => A): Vector[A] =
    given ExecutionContext = ExecutionContext.global

    Await.result(Future.traverse(ranges)((from, until) => Future(scan(from, until))), Duration.Inf)
