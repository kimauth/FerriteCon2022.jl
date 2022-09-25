abstract type Material end

#################################
# bulk material response
#################################
struct Elastic <: Material
    E::Float64
    ν::Float64
end

# for 2D only valid in plane strain
function material_response(material::Elastic, ε::SymmetricTensor{2, dim, T}) where {dim, T}
    (; E, ν) = material
    
    λ = E*ν / ((1 + ν) * (1 - 2ν))
    μ = E / (2(1 + ν))
    δ(i,j) = i == j ? one(T) : zero(T)
    f = (i,j,k,l) -> λ*δ(i,j)*δ(k,l) + μ*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k))

    C = SymmetricTensor{4, dim}(f)

    σ = C ⊡ ε
    
    return σ, C
end

#################################
# cohesive law
#################################
struct XuNeedleman <: Material
    Φₙ::Float64
    Φₜ::Float64
    δₙ::Float64
    δₜ::Float64
    Δₙˢ::Float64
end

function Φ(material::XuNeedleman, Δₙ, Δₜ::Vec)
    (; Φₙ, Φₜ, δₙ, δₜ, Δₙˢ) = material
    q = Φₜ / Φₙ
    r = Δₙˢ / δₙ
    _Φ = Φₙ + Φₙ * exp(-Δₙ/δₙ) * ((1 - r + Δₙ/δₙ) * (1-q)/(r-1) - (q + (r-q)/(r-1) * Δₙ/δₙ) * exp(-Δₜ ⋅ Δₜ / δₜ^2))
    return _Φ
end

function traction(material, Δ::Vec{dim, T}) where {dim, T}
    Δₜ = Vec{dim-1, T}(i->Δ[i])
    Δₙ = Δ[end]
    Tₙ = gradient(Δₙ -> Φ(material, Δₙ, Δₜ), Δₙ)
    Tₜ = gradient(Δₜ -> Φ(material, Δₙ, Δₜ), Δₜ)
    _T = Vec{dim, T}(i-> i<dim ? Tₜ[i] : Tₙ)
    return _T
end

function material_response(material::XuNeedleman, Δ::Vec{dim}) where dim
    dTdΔ, T = gradient(Δ -> traction(material, Δ), Δ, :all)
    return T, dTdΔ
end

