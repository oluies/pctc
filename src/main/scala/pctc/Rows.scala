package pctc

/** One row of the grandfather-paradox sweep: the CTC qubit is rotated by `RY(theta)` inside the
  * loop, and the post-selection probability has to be `cos^2(theta / 2)`.
  */
final case class GrandfatherRow(theta: Double, postSelectionProbability: Double, expected: Double):
  def absoluteError: Double = Math.abs(postSelectionProbability - expected)

/** One stage of the post-selected teleportation demonstration.
  *
  * `renormalizedProbability` is what the stage's probability becomes once the loop has renormalised
  * it, which is 1 for a loop and equal to the raw probability for a plain Bell measurement.
  */
final case class TeleportationRow(
    stage: String,
    postSelectionProbability: Double,
    renormalizedProbability: Double,
    messageFidelity: Double
)

/** One row of the retrocausal-messaging sweep. `bits` values may be infinite, which is the honest
  * answer for a noiseless CTC.
  */
final case class RetrocausalRow(
    channel: String,
    parameterName: String,
    parameter: Double,
    maxInformation: Double,
    doeblinInformation: Double,
    retrocausalCapacity: Double,
    maxInformationReference: Double,
    doeblinInformationReference: Double,
    retrocausalCapacityReference: Double,
    ordinaryClassicalCapacity: Double,
    postSelectionProbability: Double,
    evaluations: Int
):
  /** The search is accurate to about 1e-6 bits, so anything closer than that is called a tie. */
  def exceedsOrdinaryCapacity: Boolean =
    retrocausalCapacity.isInfinite || retrocausalCapacity > ordinaryClassicalCapacity + 1e-6
