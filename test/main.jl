using Ferrite, FerriteCon2022

xs = Vec.([(0., 0.), (.1, 0.), (0., .1), (.1, .1), (0., .1), (.1, .1), (0., .2), (.1, .2)])
nodes = Node.(xs)
cells = [
    Triangle((1,2,4)),
    Triangle((1,4,3)),
    CohesiveQuadrilateral((3,4,5,6)),
    Triangle((5,6,8)),
    Triangle((5,8,7)),
]
grid = Grid(cells, nodes)
addfaceset!(grid, "top", x->x[2]≈0.2)
addfaceset!(grid, "bottom", x->x[2]≈0.0)

σₘₐₓ = τₘₐₓ = 500 # MPa 
δₙ = δₜ = 0.0001 # mm 
Φₙ = σₘₐₓ * exp(1.0) * δₙ
Φₜ = sqrt(exp(1.0) / 2.0) * τₘₐₓ * δₜ
Δₙˢ = 0.1δₙ

E = 200e3 # MPa 
ν = 0.3

dim = 2
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
ip_coh_geo = MidPlaneInterpolation(ip_coh_base)
qr_coh = QuadratureRule{dim-1, RefCube}(nqp_coh)
cv_coh = CohesiveVectorValues(qr_coh, ip_coh, ip_coh_geo)

ip_bulk = Lagrange{dim, RefTetrahedron, order_f}()
ip_bulk_geo = Lagrange{dim, RefTetrahedron, order_geo}()
qr_bulk = QuadratureRule{dim, RefTetrahedron}(nqp_bulk)
cv_bulk = CellVectorValues(qr_bulk, ip_bulk, ip_bulk_geo)


dh = MixedDofHandler(grid)
field_coh = Field(:u, ip_coh, dim)
fh_coh = FieldHandler([field_coh], Set((3,)))
push!(dh, fh_coh)
field_bulk = Field(:u, ip_bulk, dim)
fh_bulk = FieldHandler([field_bulk], Set((1,2,4,5)))
push!(dh, fh_bulk)
close!(dh)

ndofs_per_cell_coh = ndofs_per_cell(dh, first(fh_coh.cellset))
nnodes_per_cell_coh = Ferrite.nnodes_per_cell(dh, first(fh_coh.cellset))
ndofs_per_cell_bulk = ndofs_per_cell(dh, first(fh_bulk.cellset))
nnodes_per_cell_bulk = Ferrite.nnodes_per_cell(dh, first(fh_bulk.cellset))

re_coh = Vector{Float64}(undef, ndofs_per_cell_coh)
ke_coh = Matrix{Float64}(undef, ndofs_per_cell_coh, ndofs_per_cell_coh)
xe_coh = Vector{Vec{dim,Float64}}(undef, nnodes_per_cell_coh)
dofs_coh = Vector{Int}(undef, ndofs_per_cell_coh)
cb_coh = CellBuffer(re_coh, ke_coh, xe_coh, dofs_coh, cv_coh)

re_bulk = Vector{Float64}(undef, ndofs_per_cell_bulk)
ke_bulk = Matrix{Float64}(undef, ndofs_per_cell_bulk, ndofs_per_cell_bulk)
xe_bulk = Vector{Vec{dim,Float64}}(undef, nnodes_per_cell_bulk)
dofs_bulk = Vector{Int}(undef, ndofs_per_cell_bulk)
cb_bulk = CellBuffer(re_coh, ke_coh, xe_bulk, dofs_bulk, cv_bulk)

subdomains = [ SubDomain(fh_coh.cellset, interface, cb_coh), SubDomain(fh_bulk.cellset, bulk, cb_bulk)]


K = create_sparsity_pattern(dh)
f = Vector{Float64}(undef, ndofs(dh))
assembler = start_assemble(K, f)

u = zeros(ndofs(dh))

uᵖ = 0.05 * 0.2 # mm 
tmax = 1.0
ch = ConstraintHandler(dh)
add!(ch, fh_bulk, Dirichlet(:u, getfaceset(grid, "top"), (x,t)->uᵖ * t/tmax, 2))
add!(ch, fh_bulk, Dirichlet(:u, getfaceset(grid, "bottom"), (x,t)->0.0, 2))
add!(ch, fh_bulk, Dirichlet(:u, Set((1,)), (x,t)->0., 1))
close!(ch)

# Newton-Raphson solver
max_iter = 20
tol = 1e-8
pvd = paraview_collection("dummy")
for (i, t) in enumerate(range(0.0, tmax; length=21))
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
                pvd[t] = vtk
            end
            break
        end
        apply_zero!(K, f, ch)
        du = K \ f
        u -= du
    end
end
vtk_save(pvd)
