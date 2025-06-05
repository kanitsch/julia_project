using Graphs
using Random
using GraphRecipes, Plots
using LinearAlgebra 


function generate_random_graph(n::Int,R_min=1.0,R_max=10.0,ε=10.0)
    g=Graph(n)
    for i in 1:n-1
        connected=rand(i+1:n)
        add_edge!(g,i,connected)
        for j in i+1:n
            if j != i && j != connected && rand() < 0.2
                add_edge!(g, i, j)
            end
        end
    end

    for i in 1:n
        while length(neighbors(g,i)) <2
            new_neighbor = rand(1:n-1)
            if new_neighbor != i && !(new_neighbor in neighbors(g, i))
                add_edge!(g, i, new_neighbor)
            end
        end                
    end

    resistances = Dict{Tuple{Int,Int}, Float64}()
    for e in edges(g)
        u, v = src(e), dst(e)
        edge_key = (min(u, v), max(u, v))
        resistances[edge_key] = rand(R_min:R_max)
    end

    all_edges = collect(keys(resistances))
    eps_edge = rand(all_edges)
    resistances[eps_edge] = ε

    return g, resistances, eps_edge, ε
end


function visualise_resistances(g)
    edgelabel_dict = Dict()
    println("Graf: ", g)
    println("Rezystancje:")
    for ((u, v), R) in resistances
        if (u,v)==(a,b)
                println("  ($u, $v): SEM $R V")
                edgelabel_dict[(u, v)]="ε=$R V"

        else    
            println("  ($u, $v): $R Ω")
            edgelabel_dict[(u, v)]="$R Ω"
        end
    end
    println("Siła elektromotoryczna ε = $Eps V między wierzchołkami $a i $b")
    graphplot(g, names=1:nv(g), edgelabel=edgelabel_dict, curves=false, nodeshape=:circle) 
end


function find_cycle_from_edge(g::Graph, a::Int, b::Int)
    visited=falses(nv(g))
    road=[a]
    stack=Int[]
    returning=false
    initial=a
    vertex=b
    while !(vertex==initial && length(road)>1)
        if !returning && vertex != initial
            visited[vertex] = true
        end
        if returning
            if isempty(stack)
                return nothing
            end
            next_v=stack[end]
            if next_v in neighbors(g,vertex)
                returning = false
                pop!(stack)
            else
                visited[vertex]=false
                pop!(road)
                next_v=road[end]
            end
        else
            returning=true
            push!(road,vertex)
            for neib in neighbors(g,vertex)
                if !visited[neib]&&(length(road)>2 || neib != initial)
                    push!(stack,neib)
                    returning = false
                end
            end

            if isempty(stack)
                return nothing
            end

            next_v=pop!(stack)

            if returning
                visited[vertex]=false
                pop!(road)
                push!(stack,next_v)
                next_v=road[end]
            end
        end
        vertex = next_v

    end
    return road
    
end

function find_n_cycles(g,n)
    cycles = []
    gcopy=deepcopy(g)
    for _ in 1:n
        edges_list = collect(edges(gcopy))
        found = false
        for i in length(edges_list):-1:1
            e=edges_list[i]
            a,b = src(e), dst(e)
            cycle = find_cycle_from_edge(gcopy,a,b)
            if cycle !== nothing
                rem_edge!(gcopy, a, b)
                push!(cycles, cycle)
                found=true
                break
            end
        end
        if !found
            break
        end
    end
    return cycles
    
end

function cycle_edges_direction(cycle,edge_list)
    result=[]
    n=length(cycle)
    for i in 1:n
        a,b=cycle[i],cycle[i == n ? 1 : i+1]
        for (idx, (u,v)) in enumerate(edge_list)
            if (u==a && v==b)
                push!(result,(idx,1))
                break
            elseif (u==b && v==a)
                push!(result,(idx,-1))
                break
            end
        end
    end
    return result    
end

function solve(g, resistances, eps_edge)
    edges_list = collect(edges(g))
    edges_list = [(src(e), dst(e)) for e in edges_list]
    e = ne(g)
    v = nv(g)
    matrix = zeros(e, e)
    vector = zeros(e)
    a, b = eps_edge

    # równania węzłowe (pomijając węzeł a - masę)
    for (k, edge) in enumerate(edges_list)
        u, w = edge
        if u != a
            matrix[u, k] = -1
        end
        if w != a
            matrix[w, k] = 1
        end
    end

    # II prawo Kirchhoffa
    cycles = find_n_cycles(g, e - v + 1)
    row=v
    for cycle in cycles
        edges_in_cycle = cycle_edges_direction(cycle, edges_list)
        if row==v
            row=a
        end
        for (idx, sign) in edges_in_cycle
            edge = edges_list[idx]
            edge_key = (min(edge...), max(edge...)) 

            if edge == eps_edge
                vector[row] += sign * resistances[edge_key]
            elseif edge == (b, a)
                vector[row] -= sign * resistances[edge_key]
            end
            matrix[row, idx] = sign * resistances[edge_key]
        end
        if row==a
            row=v+1
        else
            row+=1
        end
    end
    return matrix \ vector
end


function visualise_currents(g, sem, currents, eps_edge)
    edges_list = collect(edges(g))
    edge_currents = Dict{Tuple{Int, Int}, Float64}()
    n = nv(g)
    dg = DiGraph(n)

    edgelabels = Dict{Tuple{Int, Int}, String}()
    draw_edges = Tuple{Int, Int}[]

    for (i, e) in enumerate(edges_list)
        u, v = src(e), dst(e)
        I = currents[i]
        
        if I < 0
            I = -I
            u, v = v, u
        end

        key = (u, v)
        push!(draw_edges, key)
        edge_currents[key] = I
        add_edge!(dg, u, v)

        edgelabels[key] = "$(round(I, digits=2)) A"
        
    end

    edgelabels[eps_edge] = "SEM $(sem)V" 

    graphplot(dg;
        names=1:nv(g),
        edges=draw_edges,
        edgelabel=edgelabels,
        curves=false,
        arrow=true,
        arrowlength=7,
        arrowangle=π/8,
        nodeshape=:circle,
        title="Rozkład prądów"
    )
end

g, resistances, (a,b), Eps = generate_random_graph(8)   # a: -,    b: +
visualise_resistances(g)
currents=solve(g,resistances,(a,b))
print(currents)
visualise_currents(g,Eps, currents, (a, b))

