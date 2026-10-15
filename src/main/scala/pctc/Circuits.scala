package pctc

import scotty.quantum.density.ComplexMatrix
import scotty.quantum.math.Complex

/** Full matrices for multi-qubit gates.
  *
  * Scotty's state-vector simulator applies controlled gates without ever materialising their
  * matrices, which is the right choice there and the wrong one here: the closed-timelike-curve
  * models both need the whole unitary as an operator, to partially trace it or to build a
  * superoperator from it.
  *
  * Qubit 0 is the most significant bit, so an `n`-qubit basis index is
  * `sum over k of bit(k) << (n - 1 - k)`. That is the same ordering as
  * [[scotty.quantum.density.DensityMatrix.tensor]].
  */
object Circuits:
  val pauliX: ComplexMatrix =
    ComplexMatrix.fromRows(Seq(Seq(Complex(0), Complex(1)), Seq(Complex(1), Complex(0))))

  val hadamard: ComplexMatrix =
    val h = 1 / Math.sqrt(2)

    ComplexMatrix.fromRows(Seq(Seq(Complex(h), Complex(h)), Seq(Complex(h), Complex(-h))))

  def rotationY(theta: Double): ComplexMatrix =
    ComplexMatrix.fromRows(
      Seq(
        Seq(Complex(Math.cos(theta / 2)), Complex(-Math.sin(theta / 2))),
        Seq(Complex(Math.sin(theta / 2)), Complex(Math.cos(theta / 2)))
      )
    )

  private def bitsOf(index: Int, qubits: Int): Vector[Int] =
    Vector.tabulate(qubits)(k => (index >> (qubits - 1 - k)) & 1)

  private def indexOf(bits: Vector[Int]): Int =
    bits.indices.foldLeft(0)((acc, k) => acc | (bits(k) << (bits.length - 1 - k)))

  /** A single-qubit gate on `target`, identity elsewhere. */
  def onQubit(qubits: Int, target: Int, gate: ComplexMatrix): ComplexMatrix =
    build(qubits, (bits, column) => Seq((bits, gate, column)), target)

  /** `gate` applied to `target` when `control` is one, identity when it is zero. */
  def controlled(qubits: Int, control: Int, target: Int, gate: ComplexMatrix): ComplexMatrix =
    require(control != target, "A controlled gate needs distinct control and target qubits.")

    val dimension = 1 << qubits
    val entries = Array.fill(dimension * dimension)(Complex(0))

    for column <- 0 until dimension do
      val bits = bitsOf(column, qubits)

      if bits(control) == 0 then entries(column * dimension + column) = Complex(1)
      else
        for next <- 0 to 1 do
          val amplitude = gate(next, bits(target))

          if Complex.abs(amplitude) > 0d then
            val row = indexOf(bits.updated(target, next))

            entries(row * dimension + column) = amplitude

    ComplexMatrix(dimension, dimension, entries.toVector)

  def cnot(qubits: Int, control: Int, target: Int): ComplexMatrix =
    controlled(qubits, control, target, pauliX)

  def controlledHadamard(qubits: Int, control: Int, target: Int): ComplexMatrix =
    controlled(qubits, control, target, hadamard)

  /** Exchanges two qubits. */
  def swap(qubits: Int, first: Int, second: Int): ComplexMatrix =
    val dimension = 1 << qubits
    val entries = Array.fill(dimension * dimension)(Complex(0))

    for column <- 0 until dimension do
      val bits = bitsOf(column, qubits)
      val swapped = bits.updated(first, bits(second)).updated(second, bits(first))

      entries(indexOf(swapped) * dimension + column) = Complex(1)

    ComplexMatrix(dimension, dimension, entries.toVector)

  /** Composes gates in circuit order: the first element is applied first. */
  def inOrder(gates: ComplexMatrix*): ComplexMatrix =
    require(gates.nonEmpty, "Nothing to compose.")

    gates.reduceLeft((earlier, later) => later.multiply(earlier))

  private def build(
      qubits: Int,
      expand: (Vector[Int], Int) => Seq[(Vector[Int], ComplexMatrix, Int)],
      target: Int
  ): ComplexMatrix =
    val dimension = 1 << qubits
    val entries = Array.fill(dimension * dimension)(Complex(0))

    for column <- 0 until dimension do
      val bits = bitsOf(column, qubits)

      for (_, gate, _) <- expand(bits, column); next <- 0 to 1 do
        val amplitude = gate(next, bits(target))

        if Complex.abs(amplitude) > 0d then
          entries(indexOf(bits.updated(target, next)) * dimension + column) = amplitude

    ComplexMatrix(dimension, dimension, entries.toVector)
