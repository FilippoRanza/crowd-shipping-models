# Minimize shipment cost

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

T = zeros(stations, stations, crowd_shippers)
for i in 1:crowd_shippers
    prev = nothing
    for curr in crowd_shipper_paths[i]
        if prev !== nothing
            T[prev, curr, i] = 1
            T[curr, prev, i] = 1
        end
        prev = curr
    end
end



model = Model(Cbc.Optimizer)
set_optimizer_attribute(model, "threads", 6) 

@variable(model, x_inc[i = 1:stations, j = 1:stations, k = 1:crowd_shippers, p = 1:packages; adj_mat[i, j] == 1], Bin)
@variable(model, x_out[i = 1:stations, j = 1:stations, k = 1:crowd_shippers, p = 1:packages; adj_mat[i, j] == 1], Bin)
@variable(model, z[1:crowd_shippers], Bin)

@constraint(model, 
    start_pkg[id = 1:packages], 
    sum(T[P[id][1], j, k] .* x_inc[P[id][1], j, k, id] 
        for j = 1:stations, k = 1:crowd_shippers 
            if adj_mat[P[id][1], j] == 1)
    == 1)

@constraint(model, 
    arrive_pkg[id = 1:packages], 
    sum(T[i, P[id][2], k] .* x_out[i, P[id][2], k, id] 
        for i = 1:stations, k = 1:crowd_shippers 
            if adj_mat[i, P[id][2]] == 1) 
    == 1)


@constraint(model, 
    flow[id = 1:packages, i = 1:stations, j = 1:stations; adj_mat[i, j] == 1], 
    sum(T[i, j, k] .* x_inc[i, j, k, id] for k in 1:crowd_shippers) 
        - 
    sum(T[j, i, k] .* x_out[j, i, k, id] for k in 1:crowd_shippers) 
    == 0)


@constraint(model, 
    select_cs[k = 1:crowd_shippers], 
    sum(x_inc[i, j, k, p] 
        for i = 1:stations, j = 1:stations, p = 1:packages 
            if adj_mat[i, j] == 1) 
    <= 1000z[k])

@constraint(model, 
    one_package[i = 1:stations, j = 1:stations, k = 1:crowd_shippers; adj_mat[i, j] == 1], 
    sum(x_inc[i, j, k, p] for p in 1:packages) <= 1)

@objective(model, Min, sum(z))
optimize!(model);
if termination_status(model) == OPTIMAL
    for p in 1:packages
        println("Package $p")
        for k in 1:crowd_shippers
            println("\tCrowd Shipper $k")
            for i in 1:stations, j in 1:stations if adj_mat[i, j] == 1
                    if value(x_inc[i, j, k, p]) == 1.
                        println("\t\t$i $j");
                    end
                end
            end
        end
    end
    println(value.(z))
else    
    println("Model infeasible")
    println()
end