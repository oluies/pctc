val scala3Version = "3.9.0" // LTS
val scottyVersion = "0.7.1"
val duckdbVersion = "1.5.6.0"
val munitVersion = "1.3.6"

ThisBuild / organization := "io.github.oluies"
ThisBuild / version := "0.1.0-SNAPSHOT"
ThisBuild / scalaVersion := scala3Version

// The browser build compiles the same computation sources to JavaScript. It deliberately does not
// restructure this project into a cross-build: the only platform-specific code is `Platform`, and
// everything that cannot run in a browser -- the DuckDB layer and the command-line entry point -- is
// excluded by filename below. Scotty's sources are read from a sibling checkout, which is already a
// prerequisite for the JVM build.
lazy val jvmOnlySources = Set("Platform.scala", "Main.scala", "Results.scala", "Query.scala")

// Scotty's state-vector simulator and the parallel runner of its gate-fusion pass use parallel
// collections, which have no Scala.js artifact, and so does SemiclassicalShor, which drives both. The
// density-matrix layer the browser needs depends on none of them; the rest of the fusion pass,
// FusedKernel included, compiles for JavaScript as it is. A name that a checkout lacks matches
// nothing, so the list also covers older Scotty releases.
lazy val jvmOnlyScottySources =
  Set("QuantumSimulator.scala", "VectorWrapper.scala", "FusedRunner.scala", "SemiclassicalShor.scala")

// Scotty lives beside this repository for a local checkout and inside it in continuous
// integration, where it is checked out to ./scotty. Look in both rather than assuming one.
lazy val scottySources = Def.setting {
  val base = (ThisBuild / baseDirectory).value
  val candidates = Seq(base / "scotty", base / ".." / "scotty").map(_ / "src" / "main" / "scala")

  candidates.find(_.isDirectory).getOrElse {
    sys.error(
      "Scotty's sources were not found. The browser build compiles them directly; expected one of: " +
        candidates.map(_.getAbsolutePath).mkString(", ")
    )
  }
}

lazy val web = (project in file("web"))
  .enablePlugins(ScalaJSPlugin)
  .settings(
    name := "pctc-web",
    scalaVersion := scala3Version,
    scalacOptions ++= Seq("-deprecation", "-feature", "-unchecked", "-release:21"),
    libraryDependencies += "org.scala-js" %%% "scalajs-dom" % "2.8.1",
    Compile / unmanagedSourceDirectories ++= Seq(
      (ThisBuild / baseDirectory).value / "src" / "main" / "scala",
      scottySources.value
    ),
    Compile / unmanagedSources / excludeFilter := {
      val jvmOnly = jvmOnlySources ++ jvmOnlyScottySources
      new sbt.io.SimpleFileFilter(file =>
        jvmOnly.contains(file.getName) && !file.getPath.contains("/web/")
      )
    },
    scalaJSUseMainModuleInitializer := true,
    Compile / mainClass := Some("pctc.WebMain"),
    publish / skip := true
  )

/** Assembles the static site: the page from web/site plus the linked, optimised JavaScript. */
lazy val site = taskKey[File]("Assembles the browser demonstration into web/target/site")

lazy val root = (project in file("."))
  .settings(
    name := "pctc",
    description := "Post-selected closed timelike curve simulation on Scotty, following arXiv:2509.08965.",
    licenses := Seq("Apache-2.0" -> url("https://www.apache.org/licenses/LICENSE-2.0.txt")),
    // Scotty 0.7.1, tag v0.7.1 of oluies/scotty, installed with `sbt publishLocal publishM2` there.
    // It is also on GitHub Packages, which needs a token even to read; see Scotty's README.
    resolvers += Resolver.mavenLocal,
    libraryDependencies ++= Seq(
      "io.github.oluies" %% "scotty" % scottyVersion,
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

site := {
  val out = (ThisBuild / baseDirectory).value / "web" / "target" / "site"
  val linked = (web / Compile / fullLinkJSOutput).value
  val page = (ThisBuild / baseDirectory).value / "web" / "site"

  IO.delete(out)
  IO.createDirectory(out)
  IO.copyDirectory(page, out)
  IO.copyDirectory(linked, out)

  streams.value.log.info(s"site assembled at $out")
  out
}
