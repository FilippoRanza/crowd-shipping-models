# Minimize shipment cost

using SparseArrays;

using JuMP;
using Cbc;

include("utils.jl")

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
∞ = Inf
A = [
    1 1 ∞ 1 ∞;
    3 ∞ ∞ 3 6;
    ∞ 5 ∞ ∞ 4;
    ∞ ∞ 1 5 ∞;
    ∞ ∞ 5 ∞ 8;
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

w = ones(stations)
d = ones(stations, stations)


model = Model(Cbc.Optimizer)
set_optimizer_attribute(model, "threads", 6) 

@variable(model, x[i = 1:stations, j = 1:stations, k = 1:crowd_shippers, p = 1:packages; adj_mat[i, j] == 1 && T[k][i, j] == 1], Bin)
@variable(model, t[i = 1:stations, p = 1:packages] >= 0)
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



@expression(
    model,
    not_before_expr[k = 1:crowd_shippers, p = 1:packages, i = 1:stations],
    sum(x[i, j, k, p] for j ∈ 1:stations if adj_mat[i, j] == 1 && T[k][i, j] == 1)
)

@constraint(model,
    not_before[k = 1:crowd_shippers, p = 1:packages, i = 1:stations; not_before_expr[k, p, i] ≠ 0.0],
    t[i, p] >= A[i, k] * not_before_expr[k, p, i]
)


@expression(
    model,
    not_after_expr[p = 1:packages, i = 1:stations],
    sum(A[i, k] * x[j, i, k, p] for j ∈ 1:stations, k ∈ 1:crowd_shippers if adj_mat[j, i] == 1 && T[k][j, i] == 1)
)

@constraint(model,
    not_after[p = 1:packages, i = 1:stations; not_after_expr[p, i] ≠ 0.0 &&  P[p][2] ≠ i && P[p][1] ≠ i ],
    t[i, p] <= w[i] + not_after_expr[p, i]
)


@expression(model,
    positive_time_expr[i = 1:stations, j = 1:stations, p = 1:packages; adj_mat[i, j] == 1],
    sum(x[i, j, k, p] for k in 1:crowd_shippers if T[k][i, j] == 1)
)

@constraint(model,
    positive_time[i = 1:stations, j = 1:stations, p = 1:packages; adj_mat[i, j] == 1 && positive_time_expr[i, j, p] ≠ 0],
    t[i, p] + (w[i] + d[i, j])*positive_time_expr[i, j, p]
    <= 
    t[j, p] + 1000(1 - positive_time_expr[i, j, p])
)


@objective(model, Min, sum(z))
optimize!(model)

