package pctc

/** Renders the demonstration rows as HTML fragments.
  *
  * The command-line program formats the same rows through [[TextTable]]; a browser can do better
  * than monospace padding, so this builds real tables. Both read the same row types, so the two
  * presentations cannot drift apart in content.
  */
object WebTables:
  private def escape(text: String): String =
    text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

  private def cell(tag: String, content: String, numeric: Boolean): String =
    val attribute = if numeric then " class=\"num\"" else ""

    s"<$tag$attribute>${escape(content)}</$tag>"

  private def table(headers: Seq[String], rows: Seq[Seq[String]], numeric: Set[Int]): String =
    val head = headers.zipWithIndex.map((h, i) => cell("th", h, numeric(i))).mkString
    val body = rows
      .map(cells => cells.zipWithIndex.map((c, i) => cell("td", c, numeric(i))).mkString)
      .map(cells => s"<tr>$cells</tr>")
      .mkString

    s"<table><thead><tr>$head</tr></thead><tbody>$body</tbody></table>"

  private def number(value: Double, decimals: Int): String = TextTable.number(value, decimals)

  def grandfather(rows: Seq[GrandfatherRow]): String =
    table(
      Seq("theta / pi", "P(post-select)", "cos^2(theta/2)", "error"),
      rows.map(r =>
        Seq(
          number(r.theta / Math.PI, 4),
          number(r.postSelectionProbability, 12),
          number(r.expected, 12),
          TextTable.scientific(r.absoluteError)
        )
      ),
      Set(0, 1, 2, 3)
    )

  def teleportation(rows: Seq[TeleportationRow]): String =
    table(
      Seq("stage", "P(post-select)", "renormalised", "message fidelity"),
      rows.map(r =>
        Seq(
          r.stage,
          number(r.postSelectionProbability, 12),
          number(r.renormalizedProbability, 12),
          number(r.messageFidelity, 12)
        )
      ),
      Set(1, 2, 3)
    )

  def retrocausal(rows: Seq[RetrocausalRow]): String =
    table(
      Seq("channel", "param", "I_max", "I_doe", "C_retro", "closed form", "C_ordinary", "beats it", "P(post-select)"),
      rows.map(r =>
        Seq(
          r.channel,
          number(r.parameter, 1),
          number(r.maxInformation, 6),
          number(r.doeblinInformation, 6),
          number(r.retrocausalCapacity, 6),
          number(r.retrocausalCapacityReference, 6),
          number(r.ordinaryClassicalCapacity, 6),
          if r.exceedsOrdinaryCapacity then "yes" else "no",
          number(r.postSelectionProbability, 6)
        )
      ),
      Set(1, 2, 3, 4, 5, 6, 8)
    )

  def comparison(rows: Seq[ComparisonRow]): String =
    table(
      Seq("scenario", "Deutsch solutions", "D purity", "P(post-select)", "P purity", "T purity", "T coherent", "published in"),
      rows.map(r =>
        Seq(
          r.scenario,
          if r.deutschParameters == 0 then "unique" else s"${r.deutschParameters}-parameter family",
          number(r.deutschOutputPurity, 4),
          if r.pctcForbids then "forbidden" else number(r.pctcProbability, 4),
          if r.pctcOutputPurity.isNaN then "--" else number(r.pctcOutputPurity, 4),
          number(r.tctcOutputPurity, 4),
          number(r.tctcCoherentFraction, 4),
          r.source
        )
      ),
      Set(2, 3, 4, 5, 6)
    )
