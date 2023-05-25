using Ferrite
import FerriteViz
using FerriteCohesiveZones

# Dispatch to skip rendering of cohesive elements
FerriteViz.ntriangles(::CohesiveQuadrilateral) = 0
function FerriteViz.decompose!(coord_offset, coord_matrix, ref_coord_matrix, triangle_offset, triangle_matrix, grid, cell::CohesiveQuadrilateral)
    (coord_offset, triangle_offset)
end
Ferrite.getdim(::CohesiveQuadrilateral) = 2

# Dispatch to allow plotting cell labels for cohesive elements
FerriteViz.midpoint(cell::FerriteCohesiveZones.CohesiveCell, points) = Makie.Point2f((1/2) * (points[cell.nodes[1]] + points[cell.nodes[2]]))

# run simulation and compute displacements with corresponding dof handler
include("../main/main.jl")
solutions, dh = run_simulation("../neper/n6-id4.inp" ; order_f = 1, order_geo = 1, nqp_coh = 2, nqp_bulk = 1, tmax=2.0)

# Colors for the grains
using Colors
colors = [
    Colors.JULIA_LOGO_COLORS.purple,
    Colors.JULIA_LOGO_COLORS.red,
    Colors.JULIA_LOGO_COLORS.blue,
    Colors.JULIA_LOGO_COLORS.green,
]

# Color the grains
cellvalues = zeros(length(dh.grid.cells))
cellvalues .= NaN
cellvalues[1:16] .= 1.0
cellvalues[17:29] .= 2.0
cellvalues[30:45] .= 2.0
cellvalues[46:61] .= 3.0
cellvalues[62:75] .= 1.0
cellvalues[76:91] .= 4.0

############ Actual plotting of the animation #########
import GLMakie

# Prepare the grid for plotting
plotter = FerriteViz.MakiePlotter(dh, solutions[1])

# Fix aspect ration according to actual dimensions
factor = 5. 
nodal_umax = reshape_to_nodes(dh, factor * last(solutions), :u)
max_coordinates = [GLMakie.Point2f(node.x .+ nodal_umax[1:2, node_id]) for (node_id, node) in enumerate(dh.grid.nodes)]
min_x, max_x = extrema([c[1] for c in max_coordinates])
min_y, max_y = extrema([c[2] for c in max_coordinates])

# Prepare a figure
fig = GLMakie.Figure(resolution=(1500,1500*max_y/max_x))
ax = GLMakie.Axis(fig[1,1], limits=(min_x, max_x, min_y, max_y))
GLMakie.hidedecorations!(ax)
GLMakie.hidespines!(ax)
GLMakie.colsize!(fig.layout, 1, GLMakie.Aspect(1, max_x/max_y))
FerriteViz.cellplot!(ax, plotter,cellvalues,colormap=GLMakie.cgrad(colors, 4, categorical=true),deformation_field=:u)

# FerriteViz.wireframe!(plotter,markersize=0,strokewidth=2,deformation_field=:u,celllabels=true)
FerriteViz.wireframe!(plotter,markersize=0,strokewidth=2,deformation_field=:u)

GLMakie.record(fig,"logo.gif",[solutions[1:5:end];reverse(solutions[1:5:end])],framerate=20) do solution
    FerriteViz.update!(plotter, solution)
end
