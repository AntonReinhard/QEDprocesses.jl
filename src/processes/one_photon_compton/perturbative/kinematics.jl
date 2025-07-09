# incoming phase space layout
# alias for electron rest system (based on QEDcore.TwoBodyRestSystem)

const ComptonRestSystem{COORD} = TwoBodyRestSystem{1, COORD} where {COORD}
ComptonRestSystem(coord) = ComptonRestSystem(Val(particle_index(coord)), coord)
function ComptonRestSystem(
        run_idx::Val{1}, coord::COORD
    ) where {COORD <: QEDcore.AbstractSingleParticleCoordinate}
    throw(
        ArgumentError(
            "the first incoming particle of Compton is an electron, which has no degrees of freedom in the Compton rest frame.",
        ),
    )
end
function ComptonRestSystem(
        run_idx::Val{2}, coord::COORD
    ) where {COORD <: QEDcore.AbstractSingleParticleCoordinate{2}}
    return TwoBodyRestSystem{1}(coord)
end
ComptonRestSystem(coord::CMSEnergy) = TwoBodyRestSystem{1}(coord)
function ComptonRestSystem(::Rapidity)
    throw(
        ArgumentError(
            "there is no finite rapidity for a photon in the electron rest system"
        ),
    )
end
ComptonRestSystem() = ComptonRestSystem(Energy(2))

# outgoing phase space layout
# spherical coordinates in electron rest frame

struct ComptonSphericalLayout{INPSL <: AbstractTwoBodyInPhaseSpaceLayout} <:
    AbstractOutPhaseSpaceLayout{INPSL}
    in_psl::INPSL
end

function ComptonSphericalLayout(::TwoBodyRestSystem{2})
    throw(
        ArgumentError(
            "the first incoming particle of Compton is an electron, which has no degrees of freedom in the Compton rest frame.",
        ),
    )
end

function QEDbase.phase_space_dimension(
        proc::Compton, model::PerturbativeQED, psl::ComptonSphericalLayout
    )
    # cth, phi
    return 2
end

QEDbase.in_phase_space_layout(psl::ComptonSphericalLayout) = psl.in_psl

function QEDbase._build_momenta(
        proc::Compton,
        model::PerturbativeQED,
        psl::ComptonSphericalLayout,
        in_coords::NTuple{1, T},
        out_coords::NTuple{2, T},
    ) where {T <: Real}
    P, K = QEDbase._build_momenta(proc, model, in_phase_space_layout(psl), in_coords)
    Pt = P + K

    cth, phi = @inbounds out_coords
    sth = QEDcore.sq_diff_sqrt(one(cth), cth)
    sphi, cphi = sincos(phi)

    e = getE(Pt)
    rho_t = getMag(Pt)
    s = QEDcore.sq_diff(getE(Pt), rho_t) - 1

    x = DoubleFloat(getX(Pt))
    y = DoubleFloat(getY(Pt))
    z = DoubleFloat(getZ(Pt))
    double_t = hypot(x, y, z) * cth
    double_omega_prime = s / (2 * (e - double_t))

    Kp = double_omega_prime * SFourMomentum{T}(
        one(T), sth * cphi, sth * sphi, cth
    )
    Pp = SFourMomentum{T}(Pt - Kp)
    return (P, K), (Pp, SFourMomentum{T}(Kp))
end
