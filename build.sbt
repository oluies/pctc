val scala3Version = "3.9.0" // LTS
val scottyVersion = "0.7.0-SNAPSHOT"
val duckdbVersion = "1.5.6.0"
val munitVersion = "1.3.6"

ThisBuild / organization := "xyz.entangled"
ThisBuild / version := "0.1.0-SNAPSHOT"
ThisBuild / scalaVersion := scala3Version

lazy val root = (project in file("."))
  .settings(
    name := "pctc",
    description := "Post-selected closed timelike curve simulation on Scotty, following arXiv:2509.08965.",
    licenses := Seq("Apache-2.0" -> url("https://www.apache.org/licenses/LICENSE-2.0.txt")),
    // Scotty's `pctc` branch, installed with `sbt publishLocal publishM2` in ../scotty.
    resolvers += Resolver.mavenLocal,
    libraryDependencies ++= Seq(
      "xyz.entangled" %% "scotty" % scottyVersion,
      // The one permitted extra dependency, used only by the Main/Query layer.
      "org.duckdb" % "duckdb_jdbc" % duckdbVersion,
      "org.scalameta" %% "munit" % munitVersion % Test
    ),
    // Built and tested on the latest JDK (27 at the time of writing); `-release` keeps the emitted
    // bytecode at JDK 21. sbt accepts up to `-release:25` here, but Mill compiles on its own JVM and
    // may reject it, so both build definitions stay at 21. A native compute backend using
    // java.lang.foreign would need this raised to 25 and Mill's compiler JVM pinned to match.
    scalacOptions ++= Seq(
      "-deprecation",
      "-feature",
      "-unchecked",
      "-release:21"
    ),
    Compile / mainClass := Some("pctc.Main"),
    Compile / run / fork := true,
    Compile / run / javaOptions += "-Xss8m",
    publish / skip := true
  )
