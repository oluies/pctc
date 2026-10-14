package pctc

import scotty.quantum.density.Channel

/** A one-parameter channel family together with the closed forms it can be checked against.
  *
  * All three families are the ones arXiv:2509.08965 works out in its Supplemental Material, so every
  * numerical quantity reported for them has a reference value.
  */
final case class ChannelSpec(
    name: String,
    parameterName: String,
    build: Double => Channel,
    referenceMaxInformation: Double => Double,
    referenceDoeblinInformation: Double => Double,
    ordinaryClassicalCapacity: Double => Double,
    ordinaryCapacitySource: String
):
  def referenceRetrocausalCapacity(parameter: Double): Double =
    referenceMaxInformation(parameter) + referenceDoeblinInformation(parameter)

object ChannelSpec:
  val depolarizing: ChannelSpec = ChannelSpec(
    name = "depolarizing",
    parameterName = "p",
    build = Channel.depolarizing,
    referenceMaxInformation = ClosedForms.depolarizingMaxInformation,
    referenceDoeblinInformation = ClosedForms.depolarizingDoeblinInformation,
    ordinaryClassicalCapacity = ClosedForms.depolarizingClassicalCapacity,
    ordinaryCapacitySource = "classical capacity, King 2003 (Eq. S203)"
  )

  val amplitudeDamping: ChannelSpec = ChannelSpec(
    name = "amplitude damping",
    parameterName = "gamma",
    build = Channel.amplitudeDamping,
    referenceMaxInformation = ClosedForms.amplitudeDampingMaxInformation,
    referenceDoeblinInformation = ClosedForms.amplitudeDampingDoeblinInformation,
    ordinaryClassicalCapacity = ClosedForms.amplitudeDampingHolevoInformation(_),
    ordinaryCapacitySource = "Holevo information, Giovannetti-Fazio 2005"
  )

  val erasure: ChannelSpec = ChannelSpec(
    name = "erasure",
    parameterName = "p",
    build = Channel.erasure,
    referenceMaxInformation = ClosedForms.erasureMaxInformation,
    referenceDoeblinInformation = ClosedForms.erasureDoeblinInformation,
    ordinaryClassicalCapacity = ClosedForms.erasureClassicalCapacity,
    ordinaryCapacitySource = "classical capacity, Bennett-DiVincenzo-Smolin 1997 (Eq. S220)"
  )

  val all: Vector[ChannelSpec] = Vector(depolarizing, amplitudeDamping, erasure)
