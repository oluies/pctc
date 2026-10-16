package pctc

import scotty.quantum.density.ComplexMatrix
import scotty.quantum.math.Complex

/** A family of operators to sweep: a box of coordinate ranges plus a map from a point in that box to
  * an operator. The searches in [[Capacities]] walk the box; nothing here knows about the objective.
  */
final case class OperatorFamily(name: String, box: Vector[(Double, Double)], build: Vector[Double] => ComplexMatrix):
  def dimension: Int = box.length

/** Bloch parametrisations of the operators swept in the capacity searches.
  *
  * States use spherical coordinates on the Bloch ball, so every point of the box is a state and no
  * point is wasted. Unit-trace Hermitian operators, which the Doeblin information ranges over, use
  * the same angles with the Bloch radius mapped through `r = u / (1 - u)`. That matters: the
  * optimal `tau_B` for the amplitude-damping channel sits at Bloch radius `2 (1 + sqrt(1 - gamma)) /
  * gamma - 1`, which is 38 at `gamma = 0.1` and diverges as `gamma -> 0`. A uniform box in Bloch
  * coordinates cannot reach it; a uniform box in `u` can.
  */
object Bloch:
  /** `(t I + x X + y Y + z Z) / 2`, a Hermitian qubit operator of trace `t`. */
  def qubitOperator(t: Double, x: Double, y: Double, z: Double): ComplexMatrix =
    // Built row-major directly. `fromRows` would allocate a Seq per row and flatten them, which the
    // capacity searches would pay for a few million times over.
    ComplexMatrix(
      2,
      2,
      Vector(
        Complex((t + z) / 2),
        Complex(x / 2, -y / 2),
        Complex(x / 2, y / 2),
        Complex((t - z) / 2)
      )
    )

  /** Trace-`t` qubit operator with Bloch vector of length `radius` pointing at `(theta, phi)`. */
  def qubitSpherical(t: Double, radius: Double, theta: Double, phi: Double): ComplexMatrix =
    qubitOperator(
      t,
      radius * Math.sin(theta) * Math.cos(phi),
      radius * Math.sin(theta) * Math.sin(phi),
      radius * Math.cos(theta)
    )

  /** A unit-trace Hermitian qutrit of block-diagonal form: a qubit block of trace `t` spanning
    * `|0>, |1>`, plus `1 - t` on the erasure flag `|e>`.
    *
    * Restricting to this form is not a loss for the erasure channel. Its Choi operator commutes with
    * `I (x) (I_2 + e^{i phi} |e><e|)` for every phase, both operator inequalities are linear in the
    * swept operator, and averaging over that phase group preserves the trace as well as the
    * inequality. The optimum can therefore be taken block diagonal.
    */
  def blockQutrit(block: ComplexMatrix, blockTrace: Double): ComplexMatrix =
    val zero = Complex(0)

    ComplexMatrix(
      3,
      3,
      Vector(
        block(0, 0), block(0, 1), zero,
        block(1, 0), block(1, 1), zero,
        zero, zero, Complex(1 - blockTrace)
      )
    )

  /** `u / (1 - u)`: maps `[0, radius / (1 + radius)]` onto `[0, radius]`. */
  private def blochRadius(u: Double): Double = u / (1 - u)

  private def radiusRange(radius: Double): (Double, Double) = (0d, radius / (1 + radius))

  /** Qubit-block trace swept for the block-diagonal qutrit. Symmetric, so an odd number of grid
    * steps always includes the erasure channel's optimum at trace zero.
    */
  private val QutritTraceRange = (-2d, 2d)

  private val Angles = Vector((0d, Math.PI), (0d, 2 * Math.PI))

  /** States of the given dimension, as swept for the max-information. Every box point is a state. */
  def stateFamily(dimension: Int): OperatorFamily = dimension match
    case 2 =>
      OperatorFamily(
        "qubit Bloch ball in spherical coordinates (r, theta, phi)",
        (0d, 1d) +: Angles,
        p => qubitSpherical(1d, p(0), p(1), p(2))
      )
    case 3 =>
      OperatorFamily(
        "block-diagonal qutrit: flag weight and qubit Bloch ball (w, r, theta, phi)",
        (0d, 1d) +: (0d, 1d) +: Angles,
        p => blockQutrit(qubitSpherical(p(0), p(0) * p(1), p(2), p(3)), p(0))
      )
    case other => throw new IllegalArgumentException(s"No state family for dimension $other.")

  /** Unit-trace Hermitian operators of the given dimension, as swept for the Doeblin information.
    * These need not be positive, and the Bloch radius reaches `radius` through `u / (1 - u)`.
    */
  def hermitianFamily(dimension: Int, radius: Double): OperatorFamily = dimension match
    case 2 =>
      OperatorFamily(
        s"qubit Bloch directions with radius up to $radius (u, theta, phi)",
        radiusRange(radius) +: Angles,
        p => qubitSpherical(1d, blochRadius(p(0)), p(1), p(2))
      )
    case 3 =>
      OperatorFamily(
        s"block-diagonal qutrit: block trace in [-2, 2], Bloch radius up to $radius (t, u, theta, phi)",
        QutritTraceRange +: radiusRange(radius) +: Angles,
        p => blockQutrit(qubitSpherical(p(0), blochRadius(p(1)), p(2), p(3)), p(0))
      )
    case other => throw new IllegalArgumentException(s"No Hermitian family for dimension $other.")
