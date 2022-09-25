using GLMakie
using Colors 
include("grain_boundary_computation.jl")

# run simulation and get global variable "displacements" and "dh"
include("../main/main.jl")

#  find coordinates from grid + plot
function logo_step_data(dh, u, polygon_nodes)
    grid = dh.grid
    nodal_u = reshape_to_nodes(dh, u, :u)
    coordinates = [Point2f(node.x .+ nodal_u[1:2, node_id]) for (node_id, node) in enumerate(grid.nodes)]
    points = [[coordinates[node_idx] for node_idx in nodes_around_grain] for nodes_around_grain in polygon_nodes]
  
    mesh_coords = Vector{Point2f}[]
    for cell in grid.cells
        color = :transparent
        coords = coordinates[collect(cell.nodes)]
        if isa(cell, CohesiveQuadrilateral)
            coords .= coords[[1,2,4,3]]
        end
        push!(mesh_coords, coords)
    end

    return points, mesh_coords
end

# compute the nodes that make up the grain polygons
grid = dh.grid
cellset = grid.cellsets["COH2D4"]

min_x, max_x = extrema(node->node.x[1], grid.nodes)
min_y, max_y = extrema(node->node.x[2], grid.nodes)
addfaceset!(grid, "boundary", x-> x[1] ≈ min_x || x[1] ≈ max_x || x[2] ≈ min_y || x[2] ≈ max_y)

polygon_nodes = get_nodes_around_grains(grid, getcellset(grid, "COH2D4"), getfaceset(grid, "boundary"))

#################################
# fix aspect ration according to actual dimensions
factor = 5. 
nodal_umax = reshape_to_nodes(dh, factor * last(displacements), :u)
max_coordinates = [Point2f(node.x .+ nodal_umax[1:2, node_id]) for (node_id, node) in enumerate(grid.nodes)]
min_x, max_x = extrema([c[1] for c in max_coordinates])
min_y, max_y = extrema([c[2] for c in max_coordinates])

colors = [
    Colors.JULIA_LOGO_COLORS.purple,
    Colors.JULIA_LOGO_COLORS.red,
    Colors.JULIA_LOGO_COLORS.red,
    Colors.JULIA_LOGO_COLORS.green,
    Colors.JULIA_LOGO_COLORS.blue,
    Colors.JULIA_LOGO_COLORS.purple,
]

step = Observable(1)
u = @lift factor * displacements[$step]
plot_data = @lift logo_step_data(dh, $u, polygon_nodes)
points = @lift $plot_data[1]
mesh_coords = @lift $plot_data[2]

f = Figure(resolution=(1500,1500*max_y/max_x))
ax = Axis(f[1,1], limits=(min_x, max_x, min_y, max_y))
hidedecorations!(ax)
hidespines!(ax)

colsize!(f.layout, 1, Aspect(1, max_x/max_y))
display(f)
for i in eachindex(points[])
    line = @lift $points[i]
    poly!(ax, line; color=colors[i])
end
for i in eachindex(mesh_coords[])
    poly = @lift $mesh_coords[i]
    poly!(ax, poly; color=:transparent, strokewidth = 4)
end 

# record gif
record(f, "plotting/logo_animation.gif", eachindex(displacements); framerate = 30) do i
    step[] = i
end