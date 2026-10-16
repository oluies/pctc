package pctc

import scotty.quantum.density.{ComplexMatrix, DensityMatrix}
import scotty.quantum.math.Complex

/** One scenario evaluated under both models of a closed timelike curve.
  *
  * `deutschParameters` is the dimension of Deutsch's solution set, so zero means the self-consistent
  * state is unique and anything larger is the uniqueness ambiguity. `pctcProbability` is the weight
  * of the post-selected branch; zero means that history is forbidden outright, which is the sharpest
  * difference between the two prescriptions.
  */
final case class ComparisonRow(
    scenario: String,
    source: String,
    deutschParameters: Int,
    deutschEntropy: Double,
    deutschOutputPurity: Double,
    pctcProbability: Double,
    pctcOutputPurity: Double,
    tctcOutputPurity: Double,
    tctcCoherentFraction: Double
):
  def pctcForbids: Boolean = pctcProbability <= 1e-12

  /** Deutsch and Allen's model always admit a consistent history; post-selection sometimes does not. */
  def modelsDisagree: Boolean = pctcForbids || Math.abs(deutschOutputPurity - pctcOutputPurity) > 1e-6

  /** Allen's model adds an ordinary channel to the post-selected term, so its output is generally
    * mixed. It is not bracketed by the other two: on the state-discrimination circuit it comes out
    * more mixed than either of them.
    *
    * Allen notes the exact condition for purity below Eq. (33): the output is pure only when
    * `P|psi> = 0` or when the two terms of Eq. (33) coincide.
    */
  def tctcMixes: Boolean = tctcOutputPurity < 1d - 1e-9

/** All three models of a closed timelike curve, on the same circuits.
  *
  * The scenarios are the ones with published solutions, so the numbers in the table can be checked
  * against a paper rather than against this program's own arithmetic:
  *
  *   - the grandfather paradox, D. Deutsch, Phys. Rev. D 44, 3197 (1991), and for the post-selected
  *     side S. Lloyd et al., arXiv:1005.2219;
  *   - the unproven-theorem and state-discrimination circuits of J.-M. A. Allen, Phys. Rev. A 90,
  *     042107 (2014), arXiv:1401.4933, which works every model on each, and whose Eq. (35) gives the
  *     transition-probability output explicitly.
  */
object CtcComparison:
  private val qubit = Vector(2)
  private val root2Inverse = 1 / Math.sqrt(2)

  private def pure(dims: Vector[Int], amplitudes: Seq[Complex]): DensityMatrix =
    DensityMatrix.pure(dims, amplitudes)

  /** `Tr[rho^2]`: one for a pure state, `1/d` for the maximally mixed one. */
  def purity(state: DensityMatrix): Double = state.multiply(state).realTrace

  /** Evaluates one circuit under both models.
    *
    * `unitary` acts on the chronology-respecting system followed by the curve, in that order, and
    * the curve is the trailing `ctcDims`.
    */
  def compare(
      scenario: String,
      source: String,
      unitary: ComplexMatrix,
      crState: DensityMatrix,
      ctcDims: Vector[Int],
      grid: SearchGrid
  ): ComparisonRow =
    val deutsch = DeutschCtc.solve(unitary, crState, ctcDims, grid)

    val loopOperator = PctcLoop.loopOperator(
      unitary,
      crState.dims ++ ctcDims,
      ctcSubsystem = crState.dims.length
    )
    val transformed = loopOperator.multiply(crState.matrix).multiply(loopOperator.adjoint)
    val weight = transformed.trace.r
    val ctcDimension = ctcDims.product

    val pctcPurity =
      if weight <= 1e-12 then Double.NaN
      else purity(DensityMatrix(crState.dims, transformed.scale(1 / weight)))

    val tctc = TctcLoop.evaluate(unitary, crState, ctcDims)

    ComparisonRow(
      scenario = scenario,
      source = source,
      deutschParameters = deutsch.family.parameters,
      deutschEntropy = deutsch.entropy,
      deutschOutputPurity = purity(deutsch.output),
      pctcProbability = weight / (ctcDimension.toDouble * ctcDimension),
      pctcOutputPurity = pctcPurity,
      tctcOutputPurity = purity(tctc.normalized),
      tctcCoherentFraction = tctc.coherentFraction
    )

  /** `U = SWAP_MT CNOT_BM CNOT_TB` on `(B, M, T)`, Allen's unproven-theorem circuit. */
  val unprovenTheorem: ComplexMatrix = Circuits.inOrder(
    Circuits.cnot(qubits = 3, control = 2, target = 0),
    Circuits.cnot(qubits = 3, control = 0, target = 1),
    Circuits.swap(qubits = 3, first = 1, second = 2)
  )

  /** `U = CH_{CR,CV} SWAP`, Allen's state-discrimination circuit. */
  val discriminator: ComplexMatrix = Circuits.inOrder(
    Circuits.swap(qubits = 2, first = 0, second = 1),
    Circuits.controlledHadamard(qubits = 2, control = 0, target = 1)
  )

  /** Every scenario, in the order the table prints them. */
  def all(grid: SearchGrid = SearchGrid(steps = 13, rounds = 30, radius = 1d)): Vector[ComparisonRow] =
    val deutsch1991 = "Deutsch 1991 / Lloyd arXiv:1005.2219"
    val allen2014 = "Allen arXiv:1401.4933"

    val grandfather = Circuits.cnot(qubits = 2, control = 0, target = 1)
    val traveller = DensityMatrix.basis(qubit, 1)

    val rotations = Vector(0d, 0.5, 1d).map: fraction =>
      val theta = Math.PI * fraction

      compare(
        s"Ry(theta) on the curve, theta = ${TextTable.number(fraction, 1)} pi",
        deutsch1991,
        Circuits.onQubit(qubits = 2, target = 1, Circuits.rotationY(theta)),
        DensityMatrix.basis(qubit, 0),
        qubit,
        grid
      )

    Vector(
      compare("grandfather CNOT, traveller |1>", deutsch1991, grandfather, traveller, qubit, grid),
      compare("grandfather CNOT, traveller |0>", deutsch1991, grandfather, DensityMatrix.basis(qubit, 0), qubit, grid)
    ) ++ rotations ++ Vector(
      compare(
        "unproven theorem, |0>_B|0>_M",
        allen2014,
        unprovenTheorem,
        DensityMatrix.basis(Vector(2, 2), 0),
        qubit,
        grid
      ),
      compare("discrimination, |0>", allen2014, discriminator, DensityMatrix.basis(qubit, 0), qubit, grid),
      compare(
        "discrimination, |->",
        allen2014,
        discriminator,
        pure(qubit, Seq(Complex(root2Inverse), Complex(-root2Inverse))),
        qubit,
        grid
      ),
      compare(
        "discrimination, |+>",
        allen2014,
        discriminator,
        pure(qubit, Seq(Complex(root2Inverse), Complex(root2Inverse))),
        qubit,
        grid
      ),
      compare("discrimination, |1>", allen2014, discriminator, DensityMatrix.basis(qubit, 1), qubit, grid)
    )
