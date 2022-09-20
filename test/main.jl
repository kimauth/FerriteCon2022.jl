using Ferrite, FerriteCon2022

dim = 2
###################################
# FerriteMeshParser
###################################
using FerriteMeshParser

xymax = 0.1 # mm
node_perm = (1,2,4,3)
FerriteMeshParser.create_cell(::Type{CohesiveQuadrilateral}, node_numbers, ::FerriteMeshParser.AbaqusMeshFormat) = CohesiveQuadrilateral(ntuple(i->node_numbers[node_perm[i]], length(node_numbers)))
grid_file = "neper/n6-id4.inp"
temp_grid = get_ferrite_grid(grid_file; user_elements=Dict("COH2D4"=>CohesiveQuadrilateral))
nodes = [Node(ntuple(i->node.x[i]*xymax, dim)) for node in temp_grid.nodes]
grid = Grid(temp_grid.cells, nodes; nodesets=temp_grid.nodesets, cellsets=temp_grid.cellsets)
##################################

addfaceset!(grid, "left", x->x[1]≈0.0)
addfaceset!(grid, "bottom", x->x[2]≈0.0)
addfaceset!(grid, "top", x->x[2]≈xymax)

σₘₐₓ = τₘₐₓ = 500 # MPa 
δₙ = δₜ = 0.0005 # mm 
Φₙ = σₘₐₓ * exp(1.0) * δₙ
Φₜ = sqrt(exp(1.0) / 2.0) * τₘₐₓ * δₜ
Δₙˢ = 0.1δₙ

E = 200e3 # MPa 
ν = 0.3

order_f = 1
order_geo = 1
nqp_coh = 2
nqp_bulk = 1

thickness = 1.0

cohesive_material = XuNeedleman(Φₙ, Φₜ, δₙ, δₜ, Δₙˢ)
bulk_material = Elastic(E, ν)

interface = Interface(cohesive_material, thickness)
bulk = Bulk(bulk_material, thickness)

ip_coh_base = Lagrange{dim-1, RefCube, order_f}()
ip_coh = JumpInterpolation(ip_coh_base)
ip_coh_geo_base = Lagrange{dim-1, RefCube, order_geo}()
ip_coh_geo = MidPlaneInterpolation(ip_coh_geo_base)
qr_coh = QuadratureRule{dim-1, RefCube}(nqp_coh)
cv_coh = CohesiveVectorValues(qr_coh, ip_coh, ip_coh_geo)

ip_bulk = Lagrange{dim, RefTetrahedron, order_f}()
ip_bulk_geo = Lagrange{dim, RefTetrahedron, order_geo}()
qr_bulk = QuadratureRule{dim, RefTetrahedron}(nqp_bulk)
cv_bulk = CellVectorValues(qr_bulk, ip_bulk, ip_bulk_geo)


dh = MixedDofHandler(grid)
field_coh = Field(:u, ip_coh, dim)
fh_coh = FieldHandler([field_coh], getcellset(grid, "COH2D4"))
push!(dh, fh_coh)
field_bulk = Field(:u, ip_bulk, dim)
fh_bulk = FieldHandler([field_bulk], getcellset(grid, "CPE3"))
push!(dh, fh_bulk)
close!(dh)

ndofs_per_cell_coh = ndofs_per_cell(dh, first(fh_coh.cellset))
nnodes_per_cell_coh = Ferrite.nnodes_per_cell(dh, first(fh_coh.cellset))
ndofs_per_cell_bulk = ndofs_per_cell(dh, first(fh_bulk.cellset))
nnodes_per_cell_bulk = Ferrite.nnodes_per_cell(dh, first(fh_bulk.cellset))

re_coh = Vector{Float64}(undef, ndofs_per_cell_coh)
ke_coh = Matrix{Float64}(undef, ndofs_per_cell_coh, ndofs_per_cell_coh)
xe_coh = Vector{Ferrite.Vec{dim,Float64}}(undef, nnodes_per_cell_coh)
dofs_coh = Vector{Int}(undef, ndofs_per_cell_coh)
cb_coh = CellBuffer(re_coh, ke_coh, xe_coh, dofs_coh, cv_coh)

re_bulk = Vector{Float64}(undef, ndofs_per_cell_bulk)
ke_bulk = Matrix{Float64}(undef, ndofs_per_cell_bulk, ndofs_per_cell_bulk)
xe_bulk = Vector{Ferrite.Vec{dim,Float64}}(undef, nnodes_per_cell_bulk)
dofs_bulk = Vector{Int}(undef, ndofs_per_cell_bulk)
cb_bulk = CellBuffer(re_coh, ke_coh, xe_bulk, dofs_bulk, cv_bulk)

subdomains = [SubDomain(fh_coh.cellset, interface, cb_coh), SubDomain(fh_bulk.cellset, bulk, cb_bulk)] # , SubDomain(fh_coh_damaged.cellset, interface_weak, cb_coh)]


K = create_sparsity_pattern(dh)
f = Vector{Float64}(undef, ndofs(dh))
assembler = start_assemble(K, f)

u = zeros(ndofs(dh))

uᵖ = 0.03 * 0.1 # mm 
tmax = 1.0
ch = ConstraintHandler(dh)
add!(ch, fh_bulk, Dirichlet(:u, getfaceset(grid, "top"), (x,t)->uᵖ * t/tmax, 2))
add!(ch, fh_bulk, Dirichlet(:u, getfaceset(grid, "bottom"), (x,t)->0.0, 2))
add!(ch, fh_bulk, Dirichlet(:u, getfaceset(grid, "left"), (x,t)->0.0, 1))
close!(ch)

# Newton-Raphson solver
max_iter = 20
tol = 1e-8
pvd = paraview_collection("dummy_1el")
vm_stresses = [NaN for _ in 1:getncells(grid)]
reaction_forces = Ferrite.Vec{dim,Float64}[]
fv = FaceVectorValues(QuadratureRule{dim-1, RefTetrahedron}(2), ip_bulk, ip_bulk_geo)
nsteps = 101
for (i, t) in enumerate(range(0.0, tmax; length=nsteps))
    println("__________________________________")
    println("t=", t)
    println("__________________________________")
    update!(ch, t)
    apply!(u, ch)
    for iter = 1:max_iter
        assemble_domain!(assembler, u, subdomains, dh)
        @views residual = norm(assembler.f[ch.free_dofs])
        println("iter = ", iter, ": residual = ", residual) 
        if residual < tol
            converged = true
            vtk_grid("step_$i", grid) do vtk
                vtk_point_data(vtk, dh, u)
                vm_stresses!(vm_stresses, u, subdomains, dh)
                vtk_cell_data(vtk, vm_stresses, "vm_stress")
                pvd[t] = vtk
            end
            push!(reaction_forces, reaction_force(u, subdomains[2], dh, getfaceset(grid, "top"), fv))
            break
        end
        iter == max_iter && error("Not converged!")
        apply_zero!(K, f, ch)
        du = K \ f
        u -= du
    end
end
vtk_save(pvd)

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

lines!(range(0.0, uᵖ; length=nsteps), [rf[2] for rf in reaction_forces])


function vm_stresses!(vm_stresses, u, subdomains::AbstractVector{SubDomain}, dh)
    for subdomain in subdomains
        vm_stresses!(vm_stresses, u, subdomain, dh)
    end
    return vm_stresses
end

vm_stresses!(vm_stresses, u, subdomain::SubDomain{<:Interface}, dh) = nothing

function vm_stresses!(vm_stresses, u, subdomain::SubDomain{<:Bulk}, dh)
    (; cellset, problem, cell_buffer)  = subdomain
    (; dofs, xe, cv) = cell_buffer
    (; material) = problem
    
    for cellid in cellset
        Ferrite.cellcoords!(xe, dh, cellid)
        celldofs!(dofs, dh, cellid)
        Ferrite.reinit!(cv, xe)
        @views ue = u[dofs]
        V = 0.0
        σᵛᴹ = 0.0
        for qp in 1:getnquadpoints(cv)
            ε = function_symmetric_gradient(cv, qp, ue)
            σ, dσdε = FerriteCon2022.material_response(material, ε)
            V += getdetJdV(cv, qp)
            σᵛᴹ += sqrt(2/3*dev(σ)⊡dev(σ)) * getdetJdV(cv, qp)
        end
        vm_stresses[cellid] = σᵛᴹ / V
    end
    return vm_stresses
end