
function get_node_cell_dict(grid, cellset)
    node_cell_dict = Dict{Int, Vector{Int}}()
    for cell_idx in cellset
        cell = getcells(grid, cell_idx)
        for node in cell.nodes
            if haskey(node_cell_dict, node)
                push!(node_cell_dict[node], cell_idx)
            else
                node_cell_dict[node] = [cell_idx]
            end
        end
    end
    return node_cell_dict
end

function get_node_face_dict(grid, faceset)
    node_face_dict = Dict{Int, Vector{FaceIndex}}()
    for face_idx in faceset
        for node in getnodes(grid, face_idx)
            if haskey(node_face_dict, node)
                push!(node_face_dict[node], face_idx)
            else
                node_face_dict[node] = [face_idx]
            end
        end
    end
    return node_face_dict
end

Ferrite.getnodes(grid::Ferrite.AbstractGrid, face_idx::Ferrite.FaceIndex) = Ferrite.faces(getcells(grid, face_idx[1]))[face_idx[2]]

# for cell -> cohesive cells
function nextnode(grid, cell_idx::Int, node_idx)
    cell = getcells(grid, cell_idx)
    i = findfirst(isequal(node_idx), cell.nodes)
    if i == 1
        return cell.nodes[2]
    elseif i == 2
        return cell.nodes[1]
    elseif i == 3
        return cell.nodes[4]
    elseif i==4
        return cell.nodes[3]
    end
end

function nextnode(grid, face_idx::FaceIndex, old_node_idx)
    for node_idx in getnodes(grid, face_idx)
        if node_idx != old_node_idx
            return node_idx
        end
    end
end

function next_cell_idx(cell_idxs, old_cell_idx)
    for cell_idx in cell_idxs
        if cell_idx != old_cell_idx
            return cell_idx
        end
    end
    return nothing # no next cell available (boundary)
end

function find_nodes_around_grain!(nodes_around_grain, grid, node_cell_dict, node_idx, cell_idx, closed_loop)
    while !isnothing(cell_idx) && !closed_loop
        node_idx = nextnode(grid, cell_idx, node_idx)
        push!(nodes_around_grain, node_idx)
        if node_idx == first(nodes_around_grain)
            # closed the loop
            closed_loop = true
            break
        end
        cell_idx = next_cell_idx(node_cell_dict[node_idx], cell_idx)
    end
    return closed_loop
end

function get_nodes_for_linesegments(grid, node_cell_dict)
    nodes_around_grains = Vector{Int}[]

    while !isempty(node_cell_dict)
        nodes_around_grain = Int[]
        node_idx = first(keys(node_cell_dict)) # start node idx
        cell_idx = pop!(node_cell_dict[node_idx]) # start cell idx
        push!(nodes_around_grain, node_idx)

        closed_loop = false
        at_boundary = false

        find_nodes_around_grain!(nodes_around_grain, grid, node_cell_dict, node_idx, cell_idx, closed_loop)
        if !closed_loop && !at_boundary
            at_boundary = true
            node_idx = first(nodes_around_grain)
            cell_idx = isempty(node_cell_dict[node_idx]) ? nothing : pop!(node_cell_dict[node_idx])
            !isnothing(cell_idx) && reverse!(nodes_around_grain)
            find_nodes_around_grain!(nodes_around_grain, grid, node_cell_dict, node_idx, cell_idx, closed_loop)
        end
        push!(nodes_around_grains, nodes_around_grain)
        foreach(node_idx->delete!(node_cell_dict, node_idx), nodes_around_grain)
    end
    return nodes_around_grains
end

function get_nodes_around_grains(grid, cellset, faceset)
    node_cell_dict = get_node_cell_dict(grid, cellset)
    inner_grainboundary_nodes = get_nodes_for_linesegments(grid, node_cell_dict)

    node_face_dict = get_node_face_dict(grid, faceset)
    outer_boundary_nodes = get_nodes_for_linesegments(grid, node_face_dict)

    # sort both vectors together
    for (j, inner_node_vector) in pairs(inner_grainboundary_nodes)
        if first(inner_node_vector) != last(inner_node_vector)
            i = findfirst(nodes->first(nodes) == first(inner_node_vector) || last(nodes) == first(inner_node_vector), outer_boundary_nodes)
            if isnothing(i)
                error("No matching boundary piece for the $j-th vector.")
            end
            outer_node_vector = popat!(outer_boundary_nodes, i)
            if last(outer_node_vector) == first(inner_node_vector) && first(outer_node_vector) == last(inner_node_vector)
                append!(inner_node_vector, outer_node_vector)
            elseif last(outer_node_vector) == last(inner_node_vector) && first(outer_node_vector) == first(inner_node_vector)
                append!(inner_node_vector, reverse!(outer_node_vector))
            else
                @warn("No matching boundary piece for the $j-th vector.")
            end
        end
    end
    return inner_grainboundary_nodes
end


# ##################################
# usage example:


# grid = dh.grid
# cellset = grid.cellsets["COH2D4"]

# min_x, max_x = extrema(node->node.x[1], grid.nodes)
# min_y, max_y = extrema(node->node.x[2], grid.nodes)
# addfaceset!(grid, "boundary", x-> x[1] ≈ min_x || x[1] ≈ max_x || x[2] ≈ min_y || x[2] ≈ max_y)

# polygon_nodes = get_nodes_around_grains(grid, getcellset(grid, "COH2D4"), getfaceset(grid, "boundary"))

# #################################
# #  find coordinates from grid + plot

# points = ([[Point2(getnodes(deformed_grid, node_idx).x) for node_idx in nodes_around_grain] for nodes_around_grain in inner_grainboundary_nodes])

# f = Figure()
# ax = Axis(f[1,1])
# for line in points
#     lines!(ax, line; color=:black, linewidth = 4.)
#     # poly!(ax, line; strokewidth = 4., color=:transparent)
# end
# f
# colsize!(f.layout, 1, Aspect(1, 1))


# f = Figure()
# Axis(f[1, 1])

# polygon = Polygon(Point2f[(0, 0), (2, 0), (3, 1), (1, 1)])
# poly!(polygon, color = :transparent, strokewidth = 4)

# f