package pctc

/** A number known only to lie in `[low, high]`, refinable one bisection at a time.
  *
  * `high` is always a value the predicate accepts and `low` one it rejects, so `high` is a valid
  * upper bound at every stage and the bracket can be handed around half-finished.
  *
  * This exists so that precision is *demanded by a comparison* rather than configured ahead of one.
  * The search only ever needs to know which of two candidates is smaller, and for the overwhelming
  * majority of grid points that is settled by a single semidefinite test. Bisecting every point to a
  * fixed tolerance instead spends about forty tests on points that lose, and picking that tolerance
  * per phase is exactly where this search has had an accuracy bug before: a loose tolerance leaking
  * into the refinement rounds moved the centre, and a halving window cannot walk that back.
  */
final case class Bracket(low: Double, high: Double, accepts: Double => Boolean, tests: Int):
  def width: Double = high - low

  /** The best upper bound known so far. */
  def bound: Double = high

  /** One bisection step. */
  def refined: Bracket =
    val middle = low + width / 2

    if accepts(middle) then Bracket(low, middle, accepts, tests + 1)
    else Bracket(middle, high, accepts, tests + 1)

  /** Refines until the bracket is narrower than `tolerance`. */
  def narrowed(tolerance: Double): Bracket =
    var current = this

    while current.width > tolerance do current = current.refined

    current

object Bracket:
  /** A value already known exactly, for objectives that need no bisection at all. */
  def exact(value: Double): Bracket = Bracket(value, value, _ => true, 1)

  /** The least `x` in `[low, high]` that `accepts`, assuming `accepts` is monotone and accepts
    * `high`. Nothing is evaluated until the bracket is refined or compared.
    */
  def apply(low: Double, high: Double, accepts: Double => Boolean): Bracket =
    Bracket(low, high, accepts, 0)

  /** Which of two brackets holds the smaller value, refining only until they separate.
    *
    * Stops once one bracket lies wholly below the other, or once both are narrower than `floor`, in
    * which case the values are indistinguishable. Ties cannot be resolved by refinement, so some
    * floor is unavoidable; the point is that it is one global floor rather than a tolerance chosen
    * per phase.
    *
    * `left` wins only on a strict improvement. That is load bearing, not a detail: the spherical
    * Bloch parametrisation is degenerate at zero radius, where every angle names the same operator,
    * so a search meets many exactly tied points. Preferring the later of two ties moves the centre
    * the refinement rounds start from, and coordinate descent then walks to a different local
    * optimum -- measured as four orders of magnitude of accuracy on the amplitude-damping channel.
    */
  def lesser(left: Bracket, right: Bracket, floor: Double): (Bracket, Bracket, Boolean) =
    var a = left
    var b = right

    while a.high > b.low && b.high > a.low && (a.width > floor || b.width > floor) do
      if a.width >= b.width then a = a.refined else b = b.refined

    (a, b, a.high <= b.low || (a.width <= floor && b.width <= floor && a.high < b.high))
