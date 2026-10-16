package pctc

/** The browser has one thread available to it and no `Await`, so the coarse grid is scanned in
  * order. On the JVM this object hands the chunks to a thread pool instead; see
  * `src/main/scala/pctc/Platform.scala`.
  */
object Platform:
  val threads: Int = 1

  def mapRanges[A](ranges: Vector[(Int, Int)])(scan: (Int, Int) => A): Vector[A] =
    ranges.map((from, until) => scan(from, until))
