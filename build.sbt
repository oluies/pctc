val scala3Version = "3.3.8" // LTS
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
    scalacOptions ++= Seq(
      "-deprecation",
      "-feature",
      "-unchecked",
      "-source:3.3"
    ),
    Compile / mainClass := Some("pctc.Main"),
    Compile / run / fork := true,
    Compile / run / javaOptions += "-Xss8m",
    publish / skip := true
  )
