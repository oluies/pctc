package pctc

import java.nio.file.Files
import java.sql.Connection
import java.time.Instant

class ResultsSuite extends munit.FunSuite:
  private val grid = SearchGrid(steps = 5, rounds = 2, radius = 1e3)

  private def withTemporaryDatabase[A](body: Connection => A): A =
    val directory = Files.createTempDirectory("pctc-results")
    val path = directory.resolve("results.duckdb")

    try Results.withConnection(path.toString)(body)
    finally
      Files.deleteIfExists(path)
      Files.deleteIfExists(directory)

  private def scalar[A](connection: Connection, sql: String, read: java.sql.ResultSet => A): A =
    val statement = connection.createStatement()

    try
      val results = statement.executeQuery(sql)

      try
        results.next()
        read(results)
      finally results.close()
    finally statement.close()

  test("the schema is created idempotently and the three tables exist") {
    withTemporaryDatabase: connection =>
      Results.createSchema(connection)
      Results.createSchema(connection)

      val tables = scalar(
        connection,
        "SELECT string_agg(table_name, ',' ORDER BY table_name) FROM information_schema.tables",
        _.getString(1)
      )

      assertEquals(tables, "grandfather,retrocausal,teleportation")
  }

  test("rows round-trip, infinities included, and runs append rather than replace") {
    withTemporaryDatabase: connection =>
      Results.createSchema(connection)

      val rows = Demos.retrocausalSweep(Vector(0d, 0.5), grid)
      val first = RunContext.create("commit-one", grid)
      val second = RunContext.create("commit-two", grid.copy(steps = 7))

      assertEquals(Results.insertRetrocausal(connection, first, rows), rows.length)
      assertEquals(Results.insertRetrocausal(connection, second, rows), rows.length)

      assertEquals(scalar(connection, "SELECT count(*) FROM retrocausal", _.getLong(1)), 2L * rows.length)
      assertEquals(scalar(connection, "SELECT count(DISTINCT run_id) FROM retrocausal", _.getLong(1)), 2L)
      assertEquals(
        scalar(connection, "SELECT count(*) FROM retrocausal WHERE NOT isfinite(i_doe)", _.getLong(1)),
        6L
      )
      assertEquals(
        scalar(connection, "SELECT count(DISTINCT grid_steps) FROM retrocausal", _.getLong(1)),
        2L
      )
  }

  test("the grandfather and teleportation tables take their own rows") {
    withTemporaryDatabase: connection =>
      Results.createSchema(connection)

      val run = RunContext.create("commit", grid)

      assertEquals(Results.insertGrandfather(connection, run, Demos.grandfather(4)), 5)
      assertEquals(Results.insertTeleportation(connection, run, Demos.postselectedTeleportation), 4)

      assertEqualsDouble(
        scalar(connection, "SELECT max(absolute_error) FROM grandfather", _.getDouble(1)),
        0d,
        1e-9
      )
      assertEqualsDouble(
        scalar(connection, "SELECT min(post_selection_probability) FROM teleportation", _.getDouble(1)),
        0.25,
        1e-12
      )
  }

  test("a run context stamps a distinct id and a usable timestamp") {
    val first = RunContext.create("commit", grid)
    val second = RunContext.create("commit", grid)

    assertNotEquals(first.runId, second.runId)
    assert(!first.createdAt.isBefore(Instant.EPOCH))
    assertEquals(first.timestamp.toInstant, first.createdAt)
    assertEquals(first.grid.label, "steps=5,rounds=2,radius=1000.0")
  }
