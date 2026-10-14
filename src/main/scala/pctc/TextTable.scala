package pctc

import java.util.Locale

/** Plain-text table rendering and number formatting. Returns strings; never prints. */
object TextTable:
  /** Column alignment. Headers are aligned the same way as their column. */
  enum Align:
    case Left, Right

  def render(headers: Seq[String], aligns: Seq[Align], rows: Seq[Seq[String]]): String =
    require(aligns.length == headers.length, "One alignment per column.")
    require(rows.forall(_.length == headers.length), "Every row needs one cell per column.")

    val widths = headers.indices.map: column =>
      (headers(column) +: rows.map(_(column))).map(_.length).max

    val separator = widths.map("-" * _).mkString("-+-")

    val lines = (headers +: rows).map: row =>
      row.indices.map(column => pad(row(column), widths(column), aligns(column))).mkString(" | ")

    (lines.head +: separator +: lines.tail).mkString("\n")

  private def pad(cell: String, width: Int, align: Align): String = align match
    case Align.Left  => cell.padTo(width, ' ')
    case Align.Right => " " * (width - cell.length) + cell

  /** Fixed-point formatting that renders infinities and NaN readably. */
  def number(value: Double, decimals: Int = 6): String =
    if value.isNaN then "NaN"
    else if value.isPosInfinity then "inf"
    else if value.isNegInfinity then "-inf"
    else String.format(Locale.ROOT, s"%.${decimals}f", value)

  def scientific(value: Double, decimals: Int = 2): String =
    if value.isFinite then String.format(Locale.ROOT, s"%.${decimals}e", value) else number(value)
