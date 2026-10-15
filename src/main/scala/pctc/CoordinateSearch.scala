package pctc

import scala.concurrent.duration.Duration
import scala.concurrent.{Await, ExecutionContext, Future}

/** How finely the operator families of [[Bloch]] are swept.
  *
  * @param steps  points per axis, both in the coarse product grid and in each refinement scan
  * @param rounds coordinate-refinement rounds after the coarse grid; the window halves each round
  * @param radius half-width of the box swept for unit-trace Hermitian (not necessarily positive)
  *               operators
  */
final case class SearchGrid(steps: Int, rounds: Int, radius: Double):
  require(steps >= 2, "A grid needs at least two points per axis.")
  require(rounds >= 0, "The number of refinement rounds cannot be negative.")
  require(radius > 0d, "The Hermitian search radius has to be positive.")

  def label: String = s"steps=$steps,rounds=$rounds,radius=$radius"

object SearchGrid:
  /** Enough to reproduce the paper's closed forms to about nine decimals. */
  val default: SearchGrid = SearchGrid(steps = 9, rounds = 30, radius = 1e3)

/** Where a search ended up. `value` is `+inf` when no point in the box was feasible. */
final case class SearchResult(point: Option[Vector[Double]], value: Double, evaluations: Int)

/** A coarse product grid followed by coordinate-wise refinement in a shrinking window.
  *
  * Two things keep this cheap. The objective is given the incumbent best value and may answer `None`
  * whenever it cannot beat it, so a point that is not an improvement costs one semidefinite test
  * rather than a full bisection. And the coarse grid, which is where nearly all the points are, is
  * scanned at a loose tolerance, since all it has to do is pick a cell.
  *
  * The refinement rounds then run at full tolerance. They must: their window halves every round, so
  * a centre displaced by a loose comparison early on is a displacement the later rounds can no
  * longer walk back. Refinement points are cheap anyway, because all but a handful of them are
  * rejected by the incumbent after a single test.
  *
  * The coarse grid is scanned in parallel. Each chunk carries its own incumbent rather than sharing
  * one, which costs a little pruning and buys a deterministic result.
  */
object CoordinateSearch:
  /** `(point, incumbent, tolerance) => Some(value)` when the point beats the incumbent. */
  type Objective = (Vector[Double], Double, Double) => Option[Double]

  /** Bisection tolerance used while scanning the coarse grid, in bits. */
  val CoarseTolerance: Double = 0.125

  /** Smallest number of grid points worth handing to a thread of its own. */
  val MinimumChunk: Int = 512

  private given ExecutionContext = ExecutionContext.global

  def minimize(box: Vector[(Double, Double)], grid: SearchGrid, objective: Objective): SearchResult =
    val coarse = scanCoarseGrid(box, grid, objective)

    (0 until grid.rounds).foldLeft(coarse): (best, round) =>
      val shrink = Math.pow(0.5, round + 1)

      box.indices.foldLeft(best): (roundBest, axis) =>
        refineAxis(roundBest, box, grid.steps, shrink, OperatorOrder.LogTolerance, axis, objective)

  private def scanCoarseGrid(
      box: Vector[(Double, Double)],
      grid: SearchGrid,
      objective: Objective
  ): SearchResult =
    val perAxis = box.map((low, high) => axisValues(low, high, grid.steps))
    val radices = perAxis.map(_.length)
    val strides = radices.scanRight(1)(_ * _).tail
    val total = radices.product

    def pointAt(index: Int): Vector[Double] =
      radices.indices.map(a => perAxis(a)((index / strides(a)) % radices(a))).toVector

    def scan(from: Int, until: Int): SearchResult =
      (from until until).foldLeft(empty): (best, index) =>
        consider(best, pointAt(index), CoarseTolerance, objective)

    val threads = Runtime.getRuntime.availableProcessors
    val chunks = Math.max(1, Math.min(threads, total / MinimumChunk))

    if chunks <= 1 then scan(0, total)
    else
      val bounds = Vector.tabulate(chunks)(c => (total * c / chunks, total * (c + 1) / chunks))
      val parts = Await.result(Future.traverse(bounds)((f, u) => Future(scan(f, u))), Duration.Inf)

      parts.reduceLeft(combine)

  private def refineAxis(
      best: SearchResult,
      box: Vector[(Double, Double)],
      steps: Int,
      shrink: Double,
      tolerance: Double,
      axis: Int,
      objective: Objective
  ): SearchResult = best.point match
    case None => best
    case Some(centre) =>
      val (low, high) = box(axis)
      val half = (high - low) * shrink / 2
      val from = Math.max(low, centre(axis) - half)
      val to = Math.min(high, centre(axis) + half)

      axisValues(from, to, steps).foldLeft(best): (current, value) =>
        consider(current, centre.updated(axis, value), tolerance, objective)

  private def consider(
      best: SearchResult,
      point: Vector[Double],
      tolerance: Double,
      objective: Objective
  ): SearchResult =
    objective(point, best.value, tolerance) match
      case Some(value) if value < best.value => SearchResult(Some(point), value, best.evaluations + 1)
      case _                                 => best.copy(evaluations = best.evaluations + 1)

  /** Keeps the better of two chunk results while summing the work both of them did. */
  private def combine(left: SearchResult, right: SearchResult): SearchResult =
    val evaluations = left.evaluations + right.evaluations

    if right.value < left.value then right.copy(evaluations = evaluations)
    else left.copy(evaluations = evaluations)

  private def empty: SearchResult = SearchResult(None, OperatorOrder.Ceiling, 0)

  private def axisValues(low: Double, high: Double, steps: Int): Vector[Double] =
    if high <= low then Vector(low)
    else Vector.tabulate(steps)(i => low + (high - low) * i / (steps - 1))
