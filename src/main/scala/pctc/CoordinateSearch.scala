package pctc

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
  /** Enough to reproduce the paper's closed forms to about three decimals. */
  val default: SearchGrid = SearchGrid(steps = 9, rounds = 30, radius = 1e3)

/** Where a search ended up. `value` is `+inf` when no point in the box was feasible. */
final case class SearchResult(point: Option[Vector[Double]], value: Double, evaluations: Int)

/** A coarse product grid followed by coordinate-wise refinement in a shrinking window.
  *
  * The objective is given the incumbent best value and may answer `None` whenever it cannot beat it,
  * which lets the bisection in [[OperatorOrder]] stop after a single semidefinite test for the vast
  * majority of points.
  */
object CoordinateSearch:
  /** `(point, incumbent) => Some(value)` when the point beats the incumbent, `None` otherwise. */
  type Objective = (Vector[Double], Double) => Option[Double]

  def minimize(box: Vector[(Double, Double)], grid: SearchGrid, objective: Objective): SearchResult =
    val coarse = coarsePoints(box, grid.steps).foldLeft(SearchResult(None, OperatorOrder.Ceiling, 0)):
      (best, point) => consider(best, point, objective)

    (0 until grid.rounds).foldLeft(coarse): (best, round) =>
      val shrink = Math.pow(0.5, round + 1)

      box.indices.foldLeft(best): (roundBest, axis) =>
        refineAxis(roundBest, box, grid.steps, shrink, axis, objective)

  private def consider(best: SearchResult, point: Vector[Double], objective: Objective): SearchResult =
    objective(point, best.value) match
      case Some(value) if value < best.value => SearchResult(Some(point), value, best.evaluations + 1)
      case _                                 => best.copy(evaluations = best.evaluations + 1)

  private def refineAxis(
      best: SearchResult,
      box: Vector[(Double, Double)],
      steps: Int,
      shrink: Double,
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
        consider(current, centre.updated(axis, value), objective)

  private def coarsePoints(box: Vector[(Double, Double)], steps: Int): Iterator[Vector[Double]] =
    box
      .map((low, high) => axisValues(low, high, steps))
      .foldLeft(Iterator(Vector.empty[Double])): (acc, values) =>
        acc.flatMap(prefix => values.iterator.map(prefix :+ _))

  private def axisValues(low: Double, high: Double, steps: Int): Vector[Double] =
    if high <= low then Vector(low)
    else Vector.tabulate(steps)(i => low + (high - low) * i / (steps - 1))
