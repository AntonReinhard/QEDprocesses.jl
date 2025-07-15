using QEDcore
using Pkg
Pkg.develop(url = "/home/antonr/repos/PrecisionCarriers.jl")

using PrecisionCarriers

QEDprocesses._pert_Compton_omega_prime(rand(SFourMomentum{PrecisionCarrier{Float64}}), 0.1 * rand(PrecisionCarrier{Float64}) + 0.9)
