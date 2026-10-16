package pctc

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

  /** Fixed-point formatting that renders infinities and NaN readably.
    *
    * Goes through `java.math.BigDecimal` rather than `String.format`. `String.format` needs an
    * explicit `Locale` to avoid emitting a decimal comma under a locale like `sv_SE`, and
    * `java.util.Locale` is not implemented in Scala.js at all, so the browser build cannot link
    * against it. `BigDecimal` is implemented on both platforms and is locale-free by construction.
    *
    * It is built from `toString` rather than from the double directly, which matters: the formatter
    * rounds the shortest decimal representation of a value, not its exact binary expansion, and
    * those disagree in the last digit often enough to notice. Constructing `BigDecimal` from the
    * double itself was measurably further from the formatter than the arithmetic it replaced.
    *
    * Checked against `java.util.Formatter` over four thousand values and thirteen decimal counts.
    * The one deliberate difference: a negative value that rounds to zero prints `0.000` where the
    * formatter prints `-0.000`.
    */
  def number(value: Double, decimals: Int = 6): String =
    require(decimals >= 0, "A negative number of decimals makes no sense.")

    if value.isNaN then "NaN"
    else if value.isPosInfinity then "inf"
    else if value.isNegInfinity then "-inf"
    else
      new java.math.BigDecimal(value.toString)
        .setScale(decimals, java.math.RoundingMode.HALF_UP)
        .toPlainString

  /** Thousands-separated integer.
    *
    * Hand-rolled because `String.format`'s `%,d` reaches `java.text.DecimalFormat`, which Scala.js
    * does not implement; this module is shared with the browser build.
    */
  def grouped(value: Long): String =
    val digits = Math.abs(value).toString
    val separated = digits.reverse.grouped(3).mkString(" ").reverse

    if value < 0 then s"-$separated" else separated

  /** Scientific notation, likewise built by hand so no locale is involved. */
  def scientific(value: Double, decimals: Int = 2): String =
    if !value.isFinite then number(value)
    else if value == 0d then s"${number(0d, decimals)}e+00"
    else
      val exponent = Math.floor(Math.log10(Math.abs(value))).toInt
      val mantissa = value / Math.pow(10d, exponent.toDouble)

      // Rounding the mantissa can carry it to ten, which belongs in the exponent.
      val (shown, shifted) =
        if Math.abs(mantissa) >= 10d - Math.pow(10d, -decimals.toDouble) / 2 then (mantissa / 10d, exponent + 1)
        else (mantissa, exponent)

      val magnitude = Math.abs(shifted).toString.reverse.padTo(2, '0').reverse

      s"${number(shown, decimals)}e${if shifted < 0 then "-" else "+"}$magnitude"
