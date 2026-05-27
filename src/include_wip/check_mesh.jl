using LinearAlgebra

function check_mesh(Ωₕ::Ferrite.Grid)
    cells = Ωₕ.cells
    nodes = Ωₕ.nodes
    n_cells = length(cells)
    n_nodes = length(nodes)

    println("="^50)
    println("Grid summary")
    println("="^50)
    println("  Cells (triangles) : ", n_cells)
    println("  Nodes             : ", n_nodes)

    # -- 1. Node coordinates -----------------------------------------------
    println("\n-- Node coordinates --")
    for (i, node) in enumerate(nodes)
        println("  Node $i : ", node.x)
    end

    # -- 2. Cell connectivity ----------------------------------------------
    println("\n-- Cell connectivity --")
    for (i, cell) in enumerate(cells)
        println("  Cell $i : nodes ", cell.nodes)
    end

    # -- 3. Check all node references in cells are valid -------------------
    println("\n-- Node reference validity --")
    bad_refs = [(ci, nidx) for (ci, cell) in enumerate(cells)
                            for nidx in cell.nodes
                            if !(1 ≤ nidx ≤ n_nodes)]
    if isempty(bad_refs)
        println("  Y All cell node references are valid (1 to $n_nodes)")
    else
        println("  N Invalid references found: ", bad_refs)
    end

    # -- 4. Check for duplicate nodes -------------------------------------
    println("\n-- Duplicate node check --")
    coords = [node.x for node in nodes]
    dups = [(i,j) for i in 1:n_nodes for j in i+1:n_nodes
                  if norm(coords[i] - coords[j]) < 1e-12]
    if isempty(dups)
        println("  Y No duplicate nodes found")
    else
        println("  N Duplicate node pairs: ", dups)
    end

    # -- 5. Check for degenerate cells ------------------------------------
    println("\n-- Degenerate cell check --")
    degenerate = Int[]
    for (ci, cell) in enumerate(cells)
        q1 = nodes[cell.nodes[1]].x
        q2 = nodes[cell.nodes[2]].x
        q3 = nodes[cell.nodes[3]].x
        J  = J_Φ([q1, q2, q3])
        G  = matrix_G(J)
        area = dΩₑ(G) / 2   # triangle area = (1/2) √det(G)
        if area < 1e-14
            push!(degenerate, ci)
        end
    end
    if isempty(degenerate)
        println("  Y No degenerate cells found")
    else
        println("  N Degenerate cells (zero area): ", degenerate)
    end

    # -- 6. Cell areas and total surface area -----------------------------
    println("\n-- Cell areas --")
    total_area = 0.0
    for (ci, cell) in enumerate(cells)
        q1 = nodes[cell.nodes[1]].x
        q2 = nodes[cell.nodes[2]].x
        q3 = nodes[cell.nodes[3]].x
        J    = J_Φ([q1, q2, q3])
        G    = matrix_G(J)
        area = dΩₑ(G) / 2   # (1/2) √det(G) for a triangle
        total_area += area
        println("  Cell $ci : area = ", round(area; digits=6))
    end
    println("  Total surface area : ", round(total_area; digits=6))

    # -- 7. Euler characteristic (for closed surface: V - E + F = 2) ------
    println("\n-- Topology (Euler characteristic) --")
    edges = Set{Tuple{Int,Int}}()
    for cell in cells
        ns = sort(collect(cell.nodes))
        push!(edges, (ns[1], ns[2]))
        push!(edges, (ns[1], ns[3]))
        push!(edges, (ns[2], ns[3]))
    end
    V = n_nodes
    E = length(edges)
    F = n_cells
    χ = V - E + F
    println("  V (vertices) : $V")
    println("  E (edges)    : $E")
    println("  F (faces)    : $F")
    println("  χ = V-E+F    : $χ  ", χ == 2 ? "Y (sphere topology)" :
                                     χ == 0 ? "Y (torus topology)"  :
                                              "N (unexpected — open mesh or non-manifold?)")

    # -- 8. Node-on-surface check (user supplies expected surface) ---------
    println("\n-- Distance from expected surface --")
    println("  (Example: unit sphere — replace with your surface)")
    for (i, node) in enumerate(nodes)
        r = norm(node.x)   # distance from origin — should be 1.0 for unit sphere
        println("  Node $i : ‖x‖ = ", round(r; digits=8))
    end

    println("="^50)
end

check_mesh(Ωₕ)