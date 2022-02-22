# Minimize shipment cost

using SparseArrays;

using JuMP;
using Cbc;


P = [
    [1, 5],
    [1, 4],
    [1, 3]
]

crowd_shipper_paths = [
    [1, 2],
    [1, 3],
    [4, 5],
    [1, 2, 4],
    [3, 2, 5],
]
arrive = [
    1,
    2, 
    1, 
    2, 
    3
]

arcs = [
    (1, 2),
    (1, 3),
    (2, 3),
    (2, 4),
    (2, 5),
    (4, 5)
]

stations= 5

adj_mat = zeros(stations, stations)
for arc in arcs
    i, j = arc
    adj_mat[i, j] = 1
    adj_mat[j, i] = 1
end

packages = size(P, 1)
crowd_shippers = size(crowd_shipper_paths, 1)
arrival = ones(crowd_shippers)

T = [spzeros(Int8, stations, stations) for i in 1:crowd_shippers]
for i in 1:crowd_shippers
    prev = nothing
    for curr in crowd_shipper_paths[i]
        if prev !== nothing
            T[i][prev, curr] = 1
        end
        prev = curr
    end
end


function is_subpath(a, b)
    for i in 1:length(b)-length(a) + 1
        stat = true
        for j in 1:length(a)
            if a[j] != b[j + i - 1]
                stat = false
                break
            end
        end
        if stat
            return true
        end
    end
    false
end


time_steps = 10

model = Model(Cbc.Optimizer)
set_optimizer_attribute(model, "threads", 6) 

@variable(model, x[i = 1:stations, j = 1:stations, k = 1:crowd_shippers, p = 1:packages; adj_mat[i, j] == 1 && T[k][i, j] == 1], Bin)
@variable(model, z[1:crowd_shippers], Bin)

@constraint(model, 
    start_pkg[id = 1:packages, ], 
    sum(x[P[id][1], j, k, id] 
        for j = 1:stations, k = 1:crowd_shippers
            if adj_mat[P[id][1], j] == 1 && T[k][P[id][1], j] == 1)
    == 1)

@constraint(model, 
    arrive_pkg[id = 1:packages], 
    sum(x[j, P[id][2], k, id] 
        for j = 1:stations, k = 1:crowd_shippers
            if adj_mat[j, P[id][2]] == 1 && T[k][j, P[id][2]] == 1) 
    == 1)

@constraint(model, 
    flow[id = 1:packages, i = 1:stations; P[id][2] ≠ i && P[id][1] ≠ i], 
    sum(x[j, i, k, id] for j = 1:stations,  k = 1:crowd_shippers if adj_mat[j, i] == 1 && T[k][j, i] == 1) 
        - 
    sum(x[i, j, k, id] for j = 1:stations,  k = 1:crowd_shippers if adj_mat[i, j] == 1 && T[k][i, j] == 1) 
    == 0)

@constraint(model,
    select_cs[k = 1:crowd_shippers], 
    sum(x[i, j, k, p] 
        for i = 1:stations, j = 1:stations, p = 1:packages 
            if adj_mat[i, j] == 1 && T[k][i, j] == 1) 
    <= 1000z[k])





@constraint(model, 
    one_package[i = 1:stations, j = 1:stations, k = 1:crowd_shippers; adj_mat[i, j] == 1 && T[k][i, j] == 1], 
    sum(x[i, j, k, p] for p in 1:packages) <= 1)


@objective(model, Min, sum(z))
optimize!(model)
function print_result(model)
    if termination_status(model) == OPTIMAL
        for p in 1:packages
            println("Package $p")
            for k in 1:crowd_shippers
                output = false
        
                for i in 1:stations, j in 1:stations if adj_mat[i, j] == 1 && T[k][i, j] == 1
                    if value(x[i, j, k, p]) == 1.
                        if !output
                            println("\tCrowd Shipper $k")
                            output = true
                        end
                        println("\t\t$i -> $j");
                    end

                end
            end
        end
        end
    else    
        println("Model infeasible!")
        println()
    end         
end

