
module Usadel

using Roots

export uUsadel, uBallistic, pairbreaking

function uUsadel(Δ0, Λ, ω)
    Δd = ΔD(Λ, Δ0)
    pep = complex(-Δd^6 + 3 * Δd^4 * (Λ^2 + ω^2) + (Λ^2 + ω^2)^3 - 3 * Δd^2 * (Λ^4 - 16 * Λ^2 * ω^2 + ω^4) + 6 * Δd * Λ * sqrt(complex(-3*(Δd^2 - Λ^2)^3 * ω^2 + 9 * (Δd^4 + 7 * Δd^2 * Λ^2 + Λ^4) * ω^4 + 9 * (-Δd^2 + Λ^2) * ω^6 + 3 * ω^8)))^(1/3)
    nun = complex(ω^2 - Δd^2 + Λ^2)
    rai = sqrt(complex(ω^2 - 2 * nun / 3 + nun^2 / (3 * pep) + pep / 3))
    usa = complex(ω / (2Δd) - rai / (2Δd) + sqrt(complex(2ω^2 - 4 * nun /3 - nun^2 / (3 * pep) - pep /3 + 2 * (Δd^2 + Λ^2) * ω / rai)) / (2Δd))
    # Guarantee analyticality in the full complex plane
    usa = real(usa) + abs(imag(usa)) * 1im
    return usa
end

uBallistic(Δ0, Λ, ω) = ω/Ω(Λ, Δ0)

function Ω(Λ, Δ0)
    return (ΔD(Λ, Δ0)^(2/3) -  Λ^(2/3))^(3/2)
end

function pairbreaking(Φ, n, Δ0, ξd, R, d)
    RLP = R + d/2
    Λ = ξd^2 * Δ0 / (1.76 * π * RLP^2) * (4 * (Φ - n)^2 + d^2 / RLP^2 * (Φ^2 + (n^2)/3))
    return max(real(Λ), real(Δ0) * 1e-3)
end

function ΔD(Λ, Δ0)
    Δd = ΔΛ(real(Λ), real(Δ0))
    if !(Δd > 0)
        # Δd = ΔΛ(real(Δ0/2 - imag(ω)), real(Δ0))
        Δd = ΔΛ(real(Δ0/2 - 1e-5), real(Δ0))
    end
    return Δd
end

function P(z)
    if z <= 1
        p = π/4 * z
    else
        p = log(z + sqrt(z^2 - 1)) + z/2 * atan(1/sqrt(z^2 - 1)) - sqrt(z^2 - 1)/(2*z)
    end
    return p
end

f(Δd, Δ0 = 0.23; Λ = 0) = log(Δd/Δ0) + P(Λ/Δd)

function ΔΛ(Λ, Δ0)
    g(Δd) = f(Δd, Δ0; Λ)
    return get(find_zeros(g, 0, Δ0), 1, 0)
end

end # module Usadel
