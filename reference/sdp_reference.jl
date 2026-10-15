#!/usr/bin/env julia
"""
A second, independent reference for the max- and Doeblin informations.

`sdp_reference.py` already solves Eqs. (9) and (10) of arXiv:2509.08965 with numpy and cvxpy. This
solves the same two programs with a different modelling layer (Convex.jl) and a different solver
(SCS), so agreement between the three is evidence about the formulation rather than about one
library. The output schema matches the Python one exactly.

Both programs are linear once `lambda` is folded into the variable:

    I_max(N) = log2  min { tr Y : Y >= 0,       pi (x) Y >= N(Phi) }    with Y = lambda sigma
    I_doe(N) = -log2 max { tr W : W = W',       pi (x) W <= N(Phi) }    with W = tau / lambda

Usage:  julia --project=reference reference/sdp_reference.jl
"""

using Convex, SCS, LinearAlgebra, Printf

# SCS is a first-order solver, so its default tolerances leave about 1e-3 of slack on these
# programs. Tightening them costs milliseconds at this size and makes the fixture worth comparing
# against.
const SOLVER = Convex.MOI.OptimizerWithAttributes(
    SCS.Optimizer,
    "eps_abs" => 1e-11,
    "eps_rel" => 1e-11,
    "max_iters" => 200_000,
)

# The Doeblin program is unbounded exactly when its optimal trace is zero, and a solver reports that
# as a value a few times its own accuracy rather than as zero. The smallest legitimately finite
# trace anywhere in this sweep is 2^-8.57 = 0.0026.
const UNBOUNDED_BELOW = 1e-6

depolarizing(p) = begin
    x = ComplexF64[0 1; 1 0]
    y = ComplexF64[0 -im; im 0]
    z = ComplexF64[1 0; 0 -1]
    w = sqrt(p / 4)
    [sqrt(1 - 3p / 4) * Matrix{ComplexF64}(I, 2, 2), w * x, w * y, w * z]
end

amplitude_damping(γ) = [
    ComplexF64[1 0; 0 sqrt(1 - γ)],
    ComplexF64[0 sqrt(γ); 0 0],
]

erasure(p) = begin
    kept, lost = sqrt(1 - p), sqrt(p)
    [
        kept * ComplexF64[1 0; 0 1; 0 0],
        lost * ComplexF64[0 0; 0 0; 1 0],
        lost * ComplexF64[0 0; 0 0; 0 1],
    ]
end

const CHANNELS = [
    ("depolarizing", depolarizing),
    ("amplitude damping", amplitude_damping),
    ("erasure", erasure),
]

"N(Phi_{A'A}) on A' (x) B, with Phi the normalised maximally entangled state."
function choi(kraus)
    d_in = size(kraus[1], 2)
    φ = zeros(ComplexF64, d_in * d_in)
    for i in 0:(d_in - 1)
        φ[i * d_in + i + 1] = 1 / sqrt(d_in)
    end
    Φ = φ * φ'
    id = Matrix{ComplexF64}(I, d_in, d_in)
    sum(kron(id, k) * Φ * kron(id, k)' for k in kraus)
end

completeness_defect(kraus) =
    maximum(abs.(sum(k' * k for k in kraus) - I(size(kraus[1], 2))))

function max_information(j, d_in)
    out = size(j, 1) ÷ d_in
    π_a = Matrix{ComplexF64}(I, d_in, d_in) / d_in
    y = HermitianSemidefinite(out)
    problem = minimize(real(tr(y)), [kron(π_a, y) - j ⪰ 0])
    solve!(problem, SOLVER; silent = true)
    problem.status == Convex.MOI.OPTIMAL ? log2(problem.optval) : Inf
end

function doeblin_information(j, d_in)
    out = size(j, 1) ÷ d_in
    π_a = Matrix{ComplexF64}(I, d_in, d_in) / d_in
    w = ComplexVariable(out, out)
    problem = maximize(real(tr(w)), [w == w', j - kron(π_a, w) ⪰ 0])
    solve!(problem, SOLVER; silent = true)
    if problem.status != Convex.MOI.OPTIMAL || problem.optval <= UNBOUNDED_BELOW
        Inf
    else
        -log2(problem.optval)
    end
end

format(x) = isinf(x) ? (x > 0 ? "inf" : "-inf") : @sprintf("%.12g", x)

function main()
    out = joinpath(@__DIR__, "..", "src", "test", "resources", "reference-capacities-julia.tsv")
    rows = ["channel\tparameter\ti_max\ti_doe\tcompleteness_defect"]

    for (name, build) in CHANNELS, step in 0:10
        parameter = step / 10
        kraus = build(parameter)
        j = choi(kraus)
        d_in = size(kraus[1], 2)
        push!(rows, string(
            name, "\t", @sprintf("%.1f", parameter),
            "\t", format(max_information(j, d_in)),
            "\t", format(doeblin_information(j, d_in)),
            "\t", @sprintf("%.3e", completeness_defect(kraus)),
        ))
    end

    write(out, join(rows, "\n") * "\n")
    println("wrote $(length(rows) - 1) rows to $(normpath(out))")
    println(join(rows, "\n"))
end

main()
