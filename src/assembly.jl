struct SubDomain{P, CB}
    cellset::Set{Int}
    problem::P
    cell_buffer::CB
end

function assemble_domain!(assembler, u::AbstractVector, subdomains::AbstractVector{SubDomain}, dh::MixedDofHandler)
    fill!(assembler.K.nzval, 0.0)
    fill!(assembler.f, 0.0)
    for subdomain in subdomains
        assemble_subdomain!(assembler, u, subdomain, dh)
    end
    return nothing
end

function assemble_subdomain!(assembler, u, subdomain, dh)
    (; cellset, problem, cell_buffer) = subdomain
    (; dofs, xe, ke, re) = cell_buffer

    for cellid in cellset
        celldofs!(dofs, dh, cellid)
        Ferrite.cellcoords!(xe, dh, cellid)
        @views assemble_cell!(cell_buffer, problem, u[dofs])
        Ferrite.assemble!(assembler, dofs, ke, re)
    end
    return nothing
end