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

stations= 5
packages = size(P, 1)
crowd_shippers = size(crowd_shipper_paths, 1)

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

@variable(model, z[1:crowd_shippers], Bin)
@variable(model, x[1:stations, 1:stations, 1:crowd_shippers, 1:packages], Bin)

@constraint(model, start_pkg[id = 1:packages], sum(T[P[id][1], 1:stations, 1:crowd_shippers] .* x[P[id][1], 1:stations, 1:crowd_shippers, id])  == 1)
@constraint(model, arrive_pkg[id = 1:packages], sum(T[1:stations, P[id][2], 1:crowd_shippers] .* x[1:stations, P[id][2], 1:crowd_shippers, id])  == 1)
@constraint(model, flow[id = 1:packages, i = 1:stations, j = 1:stations], 
    sum(T[i, j, 1:crowd_shippers] .* x[i, j, 1:crowd_shippers, id]) 
        - 
    sum(T[j, i, 1:crowd_shippers] .* x[j, i, 1:crowd_shippers, id])  == 0)


@constraint(model, select_cs[k = 1:crowd_shippers], sum(x[1:stations, 1:stations, k, 1:packages]) <= 100z[k])
@constraint(model, one_package[i = 1:stations, j = 1:stations, k = 1:crowd_shippers], sum(x[i, j, k, 1:packages]) <= 1)

@objective(model, Min, sum(z))
optimize!(model);
if termination_status(model) == OPTIMAL
    for p in 1:packages
        println("Package $p")
        for k in 1:crowd_shippers
            println("\tCrowd Shipper $k")
            for i in 1:stations
                for j in 1:stations
                    if value(x[i, j, k, p]) == 1.
                        println("\t\t$i $j");
                    end
                end
            end
        end
    end
    println(value.(z))
end

