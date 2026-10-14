package pctc

import java.sql.{Connection, ResultSet}
import pctc.TextTable.Align.{Left, Right}

/** Reads `results.duckdb` back and prints the two summaries worth having across runs:
  *
  *   1. the parameter range, per channel, over which `C_retro` exceeds the ordinary classical
  *      capacity of the same channel, and
  *   2. how the searched `I_max` and `I_doe` converge on their closed forms as the grid is refined.
  *
  * Usage: `sbt "runMain pctc.Query"` or `mill pctc.runMain pctc.Query`, optionally with `--db PATH`.
  */
object Query:
  def main(args: Array[String]): Unit =
    val path = args.sliding(2).collectFirst { case Array("--db", v) => v }.getOrElse(Results.DefaultPath)

    Results.withConnection(path): connection =>
      if !hasRows(connection) then
        println(s"$path has no retrocausal rows yet. Run `sbt run` first.")
      else
        println(s"Reading $path\n")
        println(runsSection(connection))
        println(advantageSection(connection))
        println(convergenceSection(connection))

  private def hasRows(connection: Connection): Boolean =
    val tables = query(connection, "SELECT table_name FROM information_schema.tables", _.getString(1))

    tables.contains("retrocausal") &&
      query(connection, "SELECT count(*) FROM retrocausal", _.getLong(1)).headOption.exists(_ > 0)

  private def runsSection(connection: Connection): String =
    val rows = query(
      connection,
      """SELECT run_id, min(created_at)::VARCHAR, any_value(git_commit), any_value(search_resolution),
        |       count(*)
        |FROM retrocausal GROUP BY run_id ORDER BY 2""".stripMargin,
      row => Seq(row.getString(1), row.getString(2), row.getString(3), row.getString(4), row.getString(5))
    )

    s"""Runs stored in the retrocausal table
       |
       |${TextTable.render(
        Seq("run_id", "created_at", "git_commit", "search_resolution", "rows"),
        Seq(Left, Left, Left, Left, Right),
        rows
      )}
       |""".stripMargin

  private def advantageSection(connection: Connection): String =
    val rows = query(
      connection,
      """WITH latest AS (
        |  SELECT * FROM retrocausal
        |  WHERE run_id = (SELECT run_id FROM retrocausal ORDER BY created_at DESC LIMIT 1)
        |),
        |advantage AS (
        |  SELECT channel, parameter_name, parameter, c_retro, ordinary_capacity,
        |         (NOT isfinite(c_retro)) OR c_retro > ordinary_capacity + 1e-6 AS beats
        |  FROM latest
        |)
        |SELECT channel, any_value(parameter_name),
        |       min(CASE WHEN beats THEN parameter END), max(CASE WHEN beats THEN parameter END),
        |       count(*) FILTER (WHERE beats), count(*),
        |       max(CASE WHEN beats AND isfinite(c_retro) THEN c_retro - ordinary_capacity END)
        |FROM advantage GROUP BY channel ORDER BY channel""".stripMargin,
      row =>
        Seq(
          row.getString(1),
          row.getString(2),
          TextTable.number(row.getDouble(3), 2),
          TextTable.number(row.getDouble(4), 2),
          row.getString(5),
          row.getString(6),
          TextTable.number(row.getDouble(7), 4)
        )
    )

    s"""Parameter range where C_retro exceeds the ordinary classical capacity (latest run)
       |
       |${TextTable.render(
        Seq("channel", "param", "from", "to", "points", "of", "largest finite gap (bits)"),
        Seq(Left, Left, Right, Right, Right, Right, Right),
        rows
      )}
       |""".stripMargin

  private def convergenceSection(connection: Connection): String =
    val rows = query(
      connection,
      """SELECT channel, grid_steps, refinement_rounds, search_radius,
        |       count(DISTINCT run_id),
        |       max(CASE WHEN isfinite(i_max_reference) THEN abs(i_max - i_max_reference) END),
        |       max(CASE WHEN isfinite(i_doe_reference) THEN abs(i_doe - i_doe_reference) END),
        |       sum(evaluations)
        |FROM retrocausal
        |GROUP BY channel, grid_steps, refinement_rounds, search_radius
        |ORDER BY channel, grid_steps, refinement_rounds""".stripMargin,
      row =>
        Seq(
          row.getString(1),
          row.getString(2),
          row.getString(3),
          TextTable.number(row.getDouble(4), 0),
          row.getString(5),
          TextTable.scientific(row.getDouble(6)),
          TextTable.scientific(row.getDouble(7)),
          row.getString(8)
        )
    )

    s"""Convergence of the searched values on their closed forms, by grid resolution
       |    Errors are the largest absolute deviation over all finite reference values of a channel.
       |    Rows where the reference is infinite (a noiseless CTC) are excluded from the error.
       |
       |${TextTable.render(
        Seq("channel", "steps", "rounds", "radius", "runs", "max I_max error", "max I_doe error", "evaluations"),
        Seq(Left, Right, Right, Right, Right, Right, Right, Right),
        rows
      )}
       |""".stripMargin

  private def query[A](connection: Connection, sql: String, read: ResultSet => A): Vector[A] =
    val statement = connection.createStatement()

    try
      val results = statement.executeQuery(sql)

      try Iterator.continually(results).takeWhile(_.next()).map(read).toVector
      finally results.close()
    finally statement.close()
