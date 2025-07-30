include("usadel.jl")

using .Usadel
using ChunkSplitters
using Quantica
using Quantica: σ
using LinearAlgebra

default_params = (;
    ħ2ome = 76.2,
    N = 1000,
    m = 0.065,
    µ = 2,
    R = 2000)

build(; kw...) = build(merge(default_params, NamedTuple(kw)))

function build(params)
    (; R, ħ2ome, m, N, µ) = params
    a0 = R * 2sin(π / N)
    a0max = a0 * (R + a0) / R
    t = ħ2ome / (2m * a0^2)
    kF = acos(1-m*µ*a0^2/ħ2ome)/a0
    hvF = ħ2ome*sin(a0*kF)/(m*a0)
    ϵ1 = -0.5 * hvF/(2 * R)
    @show 2t, N, kF, hvF

    circle = [R * SA[cos(2π * n / N), sin(2π * n / N)] for n in 0:(N-1)]
    latN = lattice(sublat(circle), names=:N)

    SOC = @hopping((r, dr; α=5) -> -im * α * (dr[2] * σ.Ix - dr[1] * σ.zy) / (2 * a0^2); range=a0max)
    Zeeman = @onsite((; Vz=0, Vx=0, Vy=0) -> Vz * σ.zz + Vx * σ.zx + Vy * σ.Iy)

    momentum = @hopping((r, dr; n=1) -> -im*clockwise(r,dr) * peierls(r, dr, n, R)/(2a0), range = a0max)
    kineticN =
        @hopping((r, dr; n=1) -> -t * peierls(r, dr, n, R), range = a0max) +
        @onsite((; µN = µ) -> (2t - µN) * σ.zI, sublats = :N)
    Δinduced =
        @onsite((r; Δ=0.23, Δ0 = Δ, Δ1 = Δ, n0=0, n=1, ϕ=0) -> Δ0 * fluxoid(r, n0) + Δ1 * fluxoid(r, n, ϕ), sublats = :N)

    kineticN2 =
        @hopping((r, dr; n=1) -> -t * peierls2(r, dr, n, R), range = a0max) +
        @onsite((; µN = µ) -> (2t - µN) * σ.z, sublats = :N)
    Δinduced2 =
        @onsite((r; Δ=0.23, Δ0 = Δ, Δ1 = Δ, n0=0, n=1, ϕ=0) -> Δ0 * fluxoid2(r, n0) + Δ1 * fluxoid2(r, n, ϕ), sublats = :N)

    hN2 = latN |> hamiltonian(kineticN2 + Δinduced2, orbitals=2)
    hN = latN |> hamiltonian(kineticN + SOC + Δinduced + Zeeman, orbitals=4)
    p = latN |> hamiltonian(momentum, orbitals=4)


    result = (; hN2, hN, p, kF, hvF, R, t, a0, ϵ1)
    return result
end


# n is the fluxoid winding, so its twice the hopping winding: n = Φ/Φ₀ = 2πRAᵩ/(h/2e) = 2 (e/ħ) * RAᵩ => Aᵩ = n/(2R) if e = ħ = 1
Avector(n, r, R) = 0.5 * n * SA[r[2], -r[1]] / R^2
peierls(r, dr, n, R) = peierls(cis(Avector(n, r, R)' * dr))
peierls(eⁱᶲ) = SA[eⁱᶲ 0 0 0; 0 eⁱᶲ 0 0; 0 0 -conj(eⁱᶲ) 0; 0 0 0 -conj(eⁱᶲ)]
peierls2(r, dr, n, R) = peierls2(cis(Avector(n, r, R)' * dr))
peierls2(eⁱᶲ) = SA[eⁱᶲ 0; 0 -conj(eⁱᶲ)]

clockwise(r, dr) = -sign(dot(SA[-r[2], r[1]], SA[dr[1], dr[2]]))

fluxoid(r, n, ϕ=0) = σ((π/2, ϕ + π/2 + round(Int, n) * atan(r[2], r[1])), :y)
fluxoid2(r, n, ϕ=0) = σ((π/2, ϕ + round(Int, n) * atan(r[2], r[1])))

function Σfluxoid(ω, r, n, Δ0, ϕ0=0; ξ, R, d=0.0)
    nInt = round(Int, n)
    ϕ = ϕ0 + nInt * atan(r[2], r[1])
    Λ = pairbreaking(n, nInt, Δ0, ξ, R, d)
    u = uBallistic(Δ0, Λ, ω)
    Σ = -(u * σ.II - σ((π/2, ϕ), :y)) / sqrt(complex(1 - u^2))
    return Σ
end

# Σfluxoid(ω, r, n, Δ0, ϕ = 0; ξ, R, d = 0.0) = fluxoid(r, n, ϕ)

function plot_heatmap_lines(xs, ωs, ρs, b; xlabel="µN", limits=(nothing, (-0.2, 0.2)), kw...)
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel, limits, kw...)
    heatmap!(ax, xs, ωs, ρs)
    qplot!(b, color=:white)
    return fig
end

character(ψ, τ) = real(sum(x -> tr(x' * τ * x), reinterpret(SVector{4,ComplexF64}, ψ)))

## Spectrum

function sanitize_spectrum((e, psi), inds)
    sp = sortperm(e, by = real)
    e´ = view(e, sp)[inds]
    psi´ = view(psi, :, sp)[:, inds]
    for i in axes(psi´, 2)
        col = view(psi´, :, i)
        normalize!(col)
        fixphase!(col)
    end
    return e´, psi´
end

function fixphase!(ψ::AbstractVector)
    ψ ./= argmax(abs, ψ)
    return ψ
end

function pick_subspace(sp)
    ϵs, psis = sanitize_spectrum(sp[1:4, around = 0], 1:4)
    is = findall(<(0), real(ϵs))
    @assert length(is) == 2
    psi = view(psis, :, is)
    ϵ = view(ϵs, is)
    return ϵ, psi
end

## Matrix of operator on SU(2) subspace

function operator(h, o; solver = ES.ShiftInvert(ES.ArnoldiMethod(nev = 6), 0), params...)
    sp = spectrum(h; solver, params...)
    _, psi = pick_subspace(sp)
    om = o(SA[]; params...)
    return psi' * om * psi
end

## Conversion between SU(2) and SO(3)

function SU2toSO3(o::AbstractMatrix)
    σi = (0.5*tr(o * σ.I), 0.5*tr(o * σ.x), 0.5*tr(o * σ.y), 0.5*tr(o * σ.z))
    return Quantica.chopsmall.(real.(σi))
end

SU2toSO3(o::AbstractVector) = real(SA[o'*σ.x*o, o'*σ.y*o, o'*σ.z*o]) / norm(o)

function SU2toPolar(o::AbstractVector)
    v = SU2toSO3(o)
    p = hypot(v[1],v[2])
    p ≈ 0 && (p = 0.000001)
    z = sin(atan(v[3], p))
    ϕ = atan(v[2], v[1])
    return (ϕ, z)
end

## Computing (m, N, µ) from (dρ/dϵ, ħvF, kFa0, R)

function mNµ(ρ´, ħvF, kFa0; params...)
    (; R, ħ2ome) = merge(default_params, NamedTuple(params))
    m = ħ2ome*abs(cos(kFa0)/(ρ´*ħvF^3))
    a0 = abs(ρ´*ħvF^2*tan(kFa0))
    N = round(Int, π/asin(a0/2R))
    µ = -1/(ρ´*ħvF*(1+sec(kFa0)))
    return (; m, N, µ)
end

function ρ´vka(; params...)
    (; R, ħ2ome, m, N, µ) = merge(default_params, NamedTuple(params))
    a0 = R * 2sin(π / N)
    t = ħ2ome / (2m * a0^2)
    µt = µ/2t
    ρ´ = 0<=µt<=2 ? (µt-1)/((2-µt)*µt)^(3/2) * 1/(a0*t^2) : 0.0
    ħvF =  0<=µt<=2 ? 2a0*t*sqrt((2-µt)*µt) : 0.0
    kFa0 =  0<=µt<=2 ? acos(1-µt) : 0.0
    return ρ´, ħvF, kFa0
end
