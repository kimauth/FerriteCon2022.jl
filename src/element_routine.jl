# collect all data structures that need to be duplicated for parallel assembly
struct CellBuffer{dim, CV}
    re::Vector{Float64}
    ke::Matrix{Float64}
    xe::Vector{Vec{dim,Float64}}
    dofs::Vector{Int}
    cv::CV
end

abstract type Material end

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

#################################
# cohesive element routine
#################################

struct Interface{M<:Material}
    material::M
    outofdim::Float64
end

# cohesive element
function assemble_cell!(
    cell_buffer::CellBuffer,
    problem::Interface,
    ue::AbstractVector{Float64},
) 

    (; re, ke, cv, xe) = cell_buffer
    fill!(cell_buffer.re, 0.0)
    fill!(cell_buffer.ke, 0.0)
    (; material, outofdim) = problem

    # update cellvalues for current cell
    reinit!(cv, xe)

    nbase_funcs = getnbasefunctions(cv)

    for qp in 1:getnquadpoints(cv)

        # Rotation matrix
        R = get_rotation(cv, qp)
        
        dΓ = getdetJdA(cv, qp) * outofdim

        # compute separations
        Δ = function_value(cv, qp, ue)
        Δ_local = Δ ⋅ R

        T_local, dTdΔ_local = material_response(material, Δ_local)

        # rotate back to global coordinates
        T = T_local ⋅ R'
        dTdΔ = R ⋅ dTdΔ_local ⋅ R'

        # compute element force vector and element stiffness matrix
        for j in 1:nbase_funcs
            Nʲ = shape_value(cv, qp, j)
            re[j] += (T ⋅ Nʲ) * dΓ
            for i in 1:nbase_funcs
                Nⁱ = shape_value(cv, qp, i)
                ke[j,i] += dTdΔ ⋅ Nʲ ⋅ Nⁱ * dΓ
            end
        end

    end
    return nothing
end

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
# bulk element routine
#################################
struct Bulk{M<:Material}
    material::M
    outofdim::Float64
end

# continuum elements
function assemble_cell!(
    cell_buffer::CellBuffer,
    problem::Bulk,
    ue::AbstractVector{Float64},
) 

    (; re, ke, cv, xe) = cell_buffer
    fill!(cell_buffer.re, 0.0)
    fill!(cell_buffer.ke, 0.0)
    (; material, outofdim) = problem

    Ferrite.reinit!(cv, xe)

    n_basefuncs = getnbasefunctions(cv)

    for qp in 1:getnquadpoints(cv)
        ε = function_symmetric_gradient(cv, qp, ue) # strain increment
        σ, dσdε = material_response(material, ε)

        dΩ = getdetJdV(cv, qp) * outofdim
        for j in 1:n_basefuncs
            ∇Nʲ = shape_gradient(cv, qp, j)
            re[j] += ∇Nʲ ⊡ σ * dΩ 
            for i in 1:n_basefuncs
                ∇Nⁱ = shape_gradient(cv, qp, i)
                ke[i,j] += ∇Nⁱ ⊡ dσdε ⊡ ∇Nʲ * dΩ
            end
        end
    end

    return nothing
end