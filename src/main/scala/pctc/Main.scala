package pctc

import pctc.TextTable.Align.{Left, Right}
import scala.util.Try

/** Runs the three demonstrations and prints them as plain-text tables.
  *
  * Nothing here claims that information can be sent backward in time. This is a numerical model of
  * the framework of arXiv:2509.08965, and of the post-selected closed timelike curves of
  * arXiv:1005.2219 that it builds on.
  *
  * Options: `--steps N --rounds N --radius X --theta-steps N`.
  */
object Main:
  private val parameters: Vector[Double] = Vector.tabulate(11)(_ / 10d)

  def main(args: Array[String]): Unit =
    val options = Options.parse(args.toIndexedSeq)

    val grandfather = Demos.grandfather(options.thetaSteps)
    val teleportation = Demos.postselectedTeleportation
    val retrocausal = Demos.retrocausalSweep(parameters, options.grid)

    println(header(options))
    println(grandfatherTable(grandfather))
    println(teleportationTable(teleportation))
    println(retrocausalTable(retrocausal))

  private final case class Options(grid: SearchGrid, thetaSteps: Int)

  private object Options:
    def parse(args: Seq[String]): Options =
      def value(flag: String, fallback: String): String =
        args.sliding(2).collectFirst { case Seq(`flag`, v) => v }.getOrElse(fallback)

      Options(
        SearchGrid(
          steps = value("--steps", SearchGrid.default.steps.toString).toInt,
          rounds = value("--rounds", SearchGrid.default.rounds.toString).toInt,
          radius = value("--radius", SearchGrid.default.radius.toString).toDouble
        ),
        thetaSteps = value("--theta-steps", "20").toInt
      )

  /** `git rev-parse` on the working tree, so a stored row can be traced back to the code. */
  private def gitCommit: String =
    import scala.sys.process.*

    val commit = Try(Process(Seq("git", "rev-parse", "HEAD")).!!.trim).getOrElse("unknown")
    val dirty = Try(Process(Seq("git", "status", "--porcelain")).!!.trim.nonEmpty).getOrElse(false)

    if dirty then s"$commit-dirty" else commit

  private def header(options: Options): String =
    s"""P-CTC simulation on Scotty
       |arXiv:2509.08965 (Ji, Lloyd, Wilde) and arXiv:1005.2219 (Lloyd et al.)
       |
       |git_commit      $gitCommit
       |search grid     ${options.grid.label}
       |
       |This is a numerical model. No information travels backward in time.""".stripMargin

  private def grandfatherTable(rows: Seq[GrandfatherRow]): String =
    val body = rows.map: row =>
      Seq(
        TextTable.number(row.theta, 6),
        TextTable.number(row.theta / Math.PI, 4),
        TextTable.number(row.postSelectionProbability, 12),
        TextTable.number(row.expected, 12),
        TextTable.scientific(row.absoluteError)
      )

    s"""
       |(a) Grandfather paradox: RY(theta) on the CTC qubit inside a noiseless loop
       |    The post-selection probability is |Tr RY(theta) / d|^2 = cos^2(theta/2),
       |    and it vanishes at theta = pi.
       |
       |${TextTable.render(
        Seq("theta", "theta/pi", "P(post-select)", "cos^2(theta/2)", "abs error"),
        Seq(Right, Right, Right, Right, Right),
        body
      )}""".stripMargin

  private def teleportationTable(rows: Seq[TeleportationRow]): String =
    val body = rows.map: row =>
      Seq(
        row.stage,
        TextTable.number(row.postSelectionProbability, 12),
        TextTable.number(row.renormalizedProbability, 12),
        TextTable.number(row.messageFidelity, 12)
      )

    s"""
       |(b) Post-selected teleportation
       |    A Bell measurement against the future end of the pair succeeds with probability 1/4
       |    and leaves the message on the past end. Inside the loop of Eq. (1) that 1/4 is the
       |    denominator, so the renormalised probability of the desirable outcome is 1.
       |
       |${TextTable.render(
        Seq("stage", "P(post-select)", "renormalised", "message fidelity"),
        Seq(Left, Right, Right, Right),
        body
      )}""".stripMargin

  private def retrocausalTable(rows: Seq[RetrocausalRow]): String =
    val body = rows.map: row =>
      Seq(
        row.channel,
        TextTable.number(row.parameter, 1),
        TextTable.number(row.maxInformation, 6),
        TextTable.number(row.doeblinInformation, 6),
        TextTable.number(row.retrocausalCapacity, 6),
        TextTable.number(row.retrocausalCapacityReference, 6),
        TextTable.number(row.ordinaryClassicalCapacity, 6),
        if row.exceedsOrdinaryCapacity then "yes" else "no",
        TextTable.number(row.postSelectionProbability, 6)
      )

    s"""
       |(c) Retrocausal messaging through a noisy CTC
       |    I_max  = inf { log2 lambda : N(Phi) <= lambda pi (x) sigma_B }            Eq. (9)
       |    I_doe  = inf { log2 lambda : pi (x) tau_B <= lambda N(Phi) }              Eq. (10)
       |    C_retro = I_max + I_doe                                                   Eq. (11)
       |    C_ord is the ordinary classical capacity of the same channel, see ClosedForms.
       |    P(post-select) is the best a teleportation loop through this CTC can do, 2^I_max / 16.
       |
       |${TextTable.render(
        Seq(
          "channel",
          "param",
          "I_max",
          "I_doe",
          "C_retro",
          "C_retro ref",
          "C_ord",
          "beats C_ord",
          "P(post-select)"
        ),
        Seq(Left, Right, Right, Right, Right, Right, Right, Right, Right),
        body
      )}""".stripMargin
