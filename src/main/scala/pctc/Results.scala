package pctc

import java.sql.{Connection, DriverManager, PreparedStatement, Timestamp}
import java.time.Instant
import java.util.UUID

/** Everything one run of [[Main]] stamps onto its rows. */
final case class RunContext(runId: String, createdAt: Instant, gitCommit: String, grid: SearchGrid):
  def timestamp: Timestamp = Timestamp.from(createdAt)

object RunContext:
  def create(gitCommit: String, grid: SearchGrid): RunContext =
    RunContext(UUID.randomUUID().toString, Instant.now(), gitCommit, grid)

/** Append-only DuckDB storage, one table per demonstration.
  *
  * This is the only place besides [[Main]] and [[Query]] that knows a database exists; the
  * simulation itself has no idea.
  */
object Results:
  val DefaultPath: String = "results.duckdb"

  def withConnection[A](path: String)(body: Connection => A): A =
    val connection = DriverManager.getConnection(s"jdbc:duckdb:$path")

    try body(connection)
    finally connection.close()

  def createSchema(connection: Connection): Unit =
    val statement = connection.createStatement()

    try
      statement.execute(
        """CREATE TABLE IF NOT EXISTS grandfather (
          |  run_id VARCHAR, created_at TIMESTAMP, git_commit VARCHAR,
          |  demo VARCHAR, channel VARCHAR, parameter_name VARCHAR, parameter DOUBLE,
          |  post_selection_probability DOUBLE, expected_probability DOUBLE, absolute_error DOUBLE,
          |  search_resolution VARCHAR
          |)""".stripMargin
      )
      statement.execute(
        """CREATE TABLE IF NOT EXISTS teleportation (
          |  run_id VARCHAR, created_at TIMESTAMP, git_commit VARCHAR,
          |  demo VARCHAR, channel VARCHAR, stage VARCHAR,
          |  post_selection_probability DOUBLE, renormalized_probability DOUBLE,
          |  message_fidelity DOUBLE, search_resolution VARCHAR
          |)""".stripMargin
      )
      statement.execute(
        """CREATE TABLE IF NOT EXISTS retrocausal (
          |  run_id VARCHAR, created_at TIMESTAMP, git_commit VARCHAR,
          |  demo VARCHAR, channel VARCHAR, parameter_name VARCHAR, parameter DOUBLE,
          |  i_max DOUBLE, i_doe DOUBLE, c_retro DOUBLE,
          |  i_max_reference DOUBLE, i_doe_reference DOUBLE, c_retro_reference DOUBLE,
          |  ordinary_capacity DOUBLE, post_selection_probability DOUBLE,
          |  search_resolution VARCHAR, grid_steps INTEGER, refinement_rounds INTEGER,
          |  search_radius DOUBLE, evaluations INTEGER
          |)""".stripMargin
      )
    finally statement.close()

  def insertGrandfather(connection: Connection, run: RunContext, rows: Seq[GrandfatherRow]): Int =
    insert(
      connection,
      """INSERT INTO grandfather VALUES (?, ?, ?, 'grandfather paradox', 'identity (noiseless CTC)',
        |  'theta', ?, ?, ?, ?, ?)""".stripMargin,
      run,
      rows,
      (statement, row) =>
        statement.setDouble(4, row.theta)
        statement.setDouble(5, row.postSelectionProbability)
        statement.setDouble(6, row.expected)
        statement.setDouble(7, row.absoluteError)
        statement.setString(8, run.grid.label)
    )

  def insertTeleportation(connection: Connection, run: RunContext, rows: Seq[TeleportationRow]): Int =
    insert(
      connection,
      """INSERT INTO teleportation VALUES (?, ?, ?, 'post-selected teleportation',
        |  'identity (noiseless CTC)', ?, ?, ?, ?, ?)""".stripMargin,
      run,
      rows,
      (statement, row) =>
        statement.setString(4, row.stage)
        statement.setDouble(5, row.postSelectionProbability)
        statement.setDouble(6, row.renormalizedProbability)
        statement.setDouble(7, row.messageFidelity)
        statement.setString(8, run.grid.label)
    )

  def insertRetrocausal(connection: Connection, run: RunContext, rows: Seq[RetrocausalRow]): Int =
    insert(
      connection,
      """INSERT INTO retrocausal VALUES (?, ?, ?, 'retrocausal messaging', ?, ?, ?, ?, ?, ?, ?, ?, ?,
        |  ?, ?, ?, ?, ?, ?, ?)""".stripMargin,
      run,
      rows,
      (statement, row) =>
        statement.setString(4, row.channel)
        statement.setString(5, row.parameterName)
        statement.setDouble(6, row.parameter)
        statement.setDouble(7, row.maxInformation)
        statement.setDouble(8, row.doeblinInformation)
        statement.setDouble(9, row.retrocausalCapacity)
        statement.setDouble(10, row.maxInformationReference)
        statement.setDouble(11, row.doeblinInformationReference)
        statement.setDouble(12, row.retrocausalCapacityReference)
        statement.setDouble(13, row.ordinaryClassicalCapacity)
        statement.setDouble(14, row.postSelectionProbability)
        statement.setString(15, run.grid.label)
        statement.setInt(16, run.grid.steps)
        statement.setInt(17, run.grid.rounds)
        statement.setDouble(18, run.grid.radius)
        statement.setInt(19, row.evaluations)
    )

  private def insert[A](
      connection: Connection,
      sql: String,
      run: RunContext,
      rows: Seq[A],
      bind: (PreparedStatement, A) => Unit
  ): Int =
    val statement = connection.prepareStatement(sql)

    try
      rows.foreach: row =>
        statement.setString(1, run.runId)
        statement.setTimestamp(2, run.timestamp)
        statement.setString(3, run.gitCommit)
        bind(statement, row)
        statement.addBatch()

      statement.executeBatch().sum
    finally statement.close()
