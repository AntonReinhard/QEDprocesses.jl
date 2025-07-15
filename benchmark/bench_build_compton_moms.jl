using BenchmarkTools
using QEDcore
using QEDprocesses
using PrecisionCarriers
using Random


RNG = MersenneTwister(0)

PROC = Compton()
MODEL = PerturbativeQED()
PSL = ComptonSphericalLayout(ComptonRestSystem(Energy(2)))

N_PHIS = 10
N_COS = 1000
N_MOMS = N_PHIS * N_COS

PHIS_LIMS = (0.0, 2 * pi)
COS_THETAS_LIMS = (0.999, 1.0)

for FLOAT_T in [Float32, Float64]
    println("= Building $N_MOMS $FLOAT_T momenta =")

    OMEGA = FLOAT_T(1.0e3)
    PHIS = FLOAT_T.((PHIS_LIMS[2] - PHIS_LIMS[1]) .* [(0.0:(1.0 / N_PHIS):1.0)...] .+ PHIS_LIMS[1])
    COS_THETAS = FLOAT_T.((COS_THETAS_LIMS[2] - COS_THETAS_LIMS[1]) .* [(0.0:(1.0 / N_COS):1.0)...] .+ COS_THETAS_LIMS[1])

    coords = [Iterators.product(COS_THETAS, PHIS)...]
    coords = [Iterators.zip(getindex.(coords, 1), getindex.(coords, 2))...]

    bench = @benchmark moms = QEDbase._build_momenta.($PROC, $MODEL, $PSL, Ref(($OMEGA,)), $coords)
    println("Mean Time: $(trunc(mean(bench.times) / 1.0e6; digits = 3))ms")

    prec_moms = QEDbase._build_momenta.(PROC, MODEL, PSL, Ref(precify(FLOAT_T, (OMEGA,))), precify(FLOAT_T, coords))
    max_eps = 0
    for single_moms in prec_moms
        for in_out_moms in single_moms
            for pair in in_out_moms
                for mom in pair
                    max_eps = max(max_eps, PrecisionCarriers._no_epsilons.(mom)...)
                end
            end
        end
    end
    println("Maximum epsilon observed: $max_eps")

    psps = PhaseSpacePoint.(PROC, MODEL, PSL, Ref((OMEGA,)), coords)
    prec_psps = PhaseSpacePoint.(PROC, MODEL, PSL, Ref(precify(FLOAT_T, (OMEGA,))), precify(FLOAT_T, coords))

    println("")
    println("= Calculating $N_MOMS unsafe differential cross sections from coordinates, $FLOAT_T =")
    bench = @benchmark unsafe_differential_cross_section.($psps)
    println("Mean Time: $(trunc(mean(bench.times) / 1.0e6; digits = 3))ms")

    diff_cs = unsafe_differential_cross_section.(prec_psps)
    max_eps = maximum(PrecisionCarriers._no_epsilons.(diff_cs))
    println("Maximum epsilon observed: $max_eps")
    println("")
end
