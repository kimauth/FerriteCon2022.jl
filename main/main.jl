using Ferrite
using FerriteCohesiveZones
using FerriteCon2022
using FerriteMeshParser

# needed for import of cohesive elements
node_perm = (1,2,4,3)
FerriteMeshParser.create_cell(::Type{FerriteCohesiveZones.CohesiveQuadrilateral}, node_numbers, ::FerriteMeshParser.AbaqusMeshFormat) = CohesiveQuadrilateral(ntuple(i->node_numbers[node_perm[i]], length(node_numbers)))

function run_simulation(; order_f = 1, order_geo = 1, nqp_coh = 2, nqp_bulk = 1)
    
    dim = 2

    # read mesh with FerriteMeshParser
    xymax = 0.1 # mm
    grid_file = "neper/n6-id4.inp"
    temp_grid = get_ferrite_grid(grid_file; user_elements=Dict("COH2D4"=>CohesiveQuadrilateral))
    nodes = [Node(ntuple(i->node.x[i]*xymax, dim)) for node in temp_grid.nodes]
    grid = Grid(temp_grid.cells, nodes; nodesets=temp_grid.nodesets, cellsets=temp_grid.cellsets)
    # add a few facesets
    addfaceset!(grid, "left", x->x[1]≈0.0)
    addfaceset!(grid, "bottom", x->x[2]≈0.0)
    addfaceset!(grid, "top", x->x[2]≈xymax)

    # cohesive zone material
    σₘₐₓ = τₘₐₓ = 500 # MPa 
    δₙ = δₜ = 0.0005 # mm 
    Φₙ = σₘₐₓ * exp(1.0) * δₙ
    Φₜ = sqrt(exp(1.0) / 2.0) * τₘₐₓ * δₜ
    Δₙˢ = 0.1δₙ
    cohesive_material = XuNeedleman(Φₙ, Φₜ, δₙ, δₜ, Δₙˢ)

    # linear elastic material
    E = 200e3 # MPa 
    ν = 0.3
    bulk_material = Elastic(E, ν)

    thickness = 1.0

    interface = Interface(cohesive_material, thickness)
    bulk = Bulk(bulk_material, thickness)

    # cohesive Values
    ip_coh_base = Lagrange{dim-1, RefCube, order_f}()
    ip_coh = JumpInterpolation(ip_coh_base)
    ip_coh_geo_base = Lagrange{dim-1, RefCube, order_geo}()
    ip_coh_geo = MidPlaneInterpolation(ip_coh_geo_base)
    qr_coh = QuadratureRule{dim-1, RefCube}(nqp_coh)
    cv_coh = CohesiveVectorValues(qr_coh, ip_coh, ip_coh_geo)

    # bulk Values
    ip_bulk = Lagrange{dim, RefTetrahedron, order_f}()
    ip_bulk_geo = Lagrange{dim, RefTetrahedron, order_geo}()
    qr_bulk = QuadratureRule{dim, RefTetrahedron}(nqp_bulk)
    cv_bulk = CellVectorValues(qr_bulk, ip_bulk, ip_bulk_geo)

    # Dofs
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
    cb_bulk = CellBuffer(re_bulk, ke_bulk, xe_bulk, dofs_bulk, cv_bulk)

    subdomains = [SubDomain(fh_coh.cellset, interface, cb_coh), SubDomain(fh_bulk.cellset, bulk, cb_bulk)]

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
    reaction_forces = Ferrite.Vec{dim,Float64}[]
    fv = FaceVectorValues(QuadratureRule{dim-1, RefTetrahedron}(2), ip_bulk, ip_bulk_geo)
    nsteps = 101
    displacements = Vector{Float64}[]
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
                push!(reaction_forces, reaction_force(u, subdomains[2], dh, getfaceset(grid, "top"), fv))
                push!(displacements, copy(u))
                break
            end
            iter == max_iter && error("Not converged!")
            apply_zero!(K, f, ch)
            du = K \ f
            u -= du

        end
    end
    return displacements, dh
end

displacements, dh = run_simulation()
