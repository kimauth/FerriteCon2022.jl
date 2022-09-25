# collect all data structures that need to be duplicated for parallel assembly
struct CellBuffer{dim, CV}
    re::Vector{Float64}
    ke::Matrix{Float64}
    xe::Vector{Vec{dim,Float64}}
    dofs::Vector{Int}
    cv::CV
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

# reaction force computation 
function reaction_force(u, subdomain, dh, faceset, fv::FaceValues{dim}) where dim
    (; problem, cell_buffer)  = subdomain
    (; dofs, xe) = cell_buffer
    (; material) = problem
    
    rf = zero(Ferrite.Vec{dim})
    for (cellid, faceid) in faceset
        Ferrite.cellcoords!(xe, dh, cellid)
        celldofs!(dofs, dh, cellid)
        Ferrite.reinit!(fv, xe, faceid)
        @views ue = u[dofs]
        for qp in 1:getnquadpoints(fv)
            n = getnormal(fv, qp)
            ε = function_symmetric_gradient(fv, qp, ue)
            σ, dσdε = FerriteCon2022.material_response(material, ε)
            dV = getdetJdV(fv, qp)
            rf += σ ⋅ n * dV
        end
    end
    return rf 
end
