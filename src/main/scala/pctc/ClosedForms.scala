package pctc

/** Closed-form reference values, used to check the numerical searches and to put the retrocausal
  * capacity next to the ordinary classical capacity of the same channel.
  *
  * Equation numbers in `Sxxx` form refer to the Supplemental Material of arXiv:2509.08965v3,
  * "Retrocausal capacity of a quantum channel", Ji, Lloyd and Wilde. All of these are for a qubit
  * input, `d_A = 2`.
  *
  * The six max- and Doeblin informations are proved, as the optimal values of Eqs. (9) and (10)
  * with `lambda` folded into the variable, in `lean/Pctc/ClosedForms.lean`:
  * `Pctc.depolarizing_maxInformation` and its five siblings, named after the functions here.
  */
object ClosedForms:
  private val Log2 = Math.log(2)

  def log2(x: Double): Double = Math.log(x) / Log2

  /** `h_2(x) = -x log2 x - (1 - x) log2 (1 - x)`. */
  def binaryEntropy(x: Double): Double =
    if x <= 0d || x >= 1d then 0d else -x * log2(x) - (1 - x) * log2(1 - x)

  // -- Max-information ------------------------------------------------------------------------

  /** `I_max(D_p) = log2[(1 - p) d^2 + p]`, Eq. (S199). */
  def depolarizingMaxInformation(p: Double): Double = log2(4 * (1 - p) + p)

  /** `I_max(E_p) = log2[(1 - p) d^2 + p]`, Eq. (S216). */
  def erasureMaxInformation(p: Double): Double = log2(4 * (1 - p) + p)

  /** `I_max(A_gamma) = log2[2 (1 + sqrt(1 - gamma)) - gamma]`, Eq. (S223). */
  def amplitudeDampingMaxInformation(gamma: Double): Double =
    log2(2 * (1 + Math.sqrt(1 - gamma)) - gamma)

  // -- Doeblin information --------------------------------------------------------------------

  /** `I_doe(D_p) = I_doe^inf(D_p) = -log2 p`, Eq. (S198). */
  def depolarizingDoeblinInformation(p: Double): Double = -log2(p)

  /** `I_doe(E_p) = I_doe^inf(E_p) = -log2 p`, Eq. (S215). */
  def erasureDoeblinInformation(p: Double): Double = -log2(p)

  /** `I_doe(A_gamma) = -2 log2[1 - sqrt(1 - gamma)]`, Eq. (S222). Only a lower bound on the
    * regularised Doeblin information, which is not known in closed form for this channel.
    */
  def amplitudeDampingDoeblinInformation(gamma: Double): Double =
    -2 * log2(1 - Math.sqrt(1 - gamma))

  // -- Ordinary (unassisted) classical capacities ---------------------------------------------

  /** `C(D_p) = log2 d - h_2(p (d - 1) / d) - (p (d - 1) / d) log2 (d - 1)`, Eq. (S203), which for a
    * qubit is `1 - h_2(p / 2)`. Proven additive for the depolarizing channel by C. King, IEEE Trans.
    * Inf. Theory 49, 221 (2003).
    */
  def depolarizingClassicalCapacity(p: Double): Double = 1 - binaryEntropy(p / 2)

  /** `C(E_p) = (1 - p) log2 d`, Eq. (S220); for a qubit, `1 - p`. C. H. Bennett, D. P. DiVincenzo and
    * J. A. Smolin, Phys. Rev. Lett. 78, 3217 (1997).
    */
  def erasureClassicalCapacity(p: Double): Double = 1 - p

  /** Holevo information of the amplitude-damping channel,
    * `chi(A_gamma) = max_{t in [0, 1]} [h_2((1 - gamma) t) - t h_2(gamma)]`, attained on the
    * ensemble `{(1 - t, |0>), (t, |1>)}`. V. Giovannetti and R. Fazio, Phys. Rev. A 71, 032314
    * (2005).
    *
    * This is the single-letter quantity. The paper quotes quantum, not classical, capacities for
    * this channel, so it is reported here as the standard reference value for comparison and is in
    * any case a lower bound on the classical capacity.
    */
  def amplitudeDampingHolevoInformation(gamma: Double, steps: Int = 100000): Double =
    val noise = binaryEntropy(gamma)

    (0 to steps).foldLeft(0d): (best, i) =>
      val t = i.toDouble / steps

      Math.max(best, binaryEntropy((1 - gamma) * t) - t * noise)
