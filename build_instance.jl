struct ModelConfig
    stations::Int
    crowd_shippers::Int
    packages::Int
end

struct NetworkInfo
    adj::Matrix{Int8}
    ct::Vector{SparseMatrixCSC{Int8,Int64}}
    wait::Vector{Float64}
    delay::Matrix{Float64}
end


function build_model(model::Model, conf::ModelConfig, net::NetworkInfo, parcels, arrival)
    @variable(
        model,
        x[
            i = 1:conf.stations,
            j = 1:conf.stations,
            k = 1:conf.crowd_shippers,
            p = 1:conf.packages;
            csr_direct_path(net, i, j, k),
        ],
        Bin
    )

    @variable(model, t[i = 1:conf.stations, p = 1:conf.packages] >= 0)
    @variable(model, z[1:conf.crowd_shippers], Bin)



    @constraint(
        model,
        start_pkg[id = 1:conf.packages,],
        sum(
            x[parcels[id, 1], j, k, id] for
            j = 1:conf.stations, k = 1:conf.crowd_shippers if csr_direct_path(net, parcels[id, 1], j, k)
        ) == 1
    )

    @constraint(
        model,
        arrive_pkg[id = 1:conf.packages],
        sum(
            x[j, parcels[id, 2], k, id] for
            j = 1:conf.stations, k = 1:conf.crowd_shippers if csr_direct_path(net, j, parcels[id, 2], k)
        ) == 1
    )


    @constraint(
        model,
        flow[id = 1:conf.packages, i = 1:conf.stations; parcels[id, 2] ≠ i && parcels[id, 1] ≠ i],
        sum(
            x[j, i, k, id] for
            j = 1:conf.stations, k = 1:conf.crowd_shippers if csr_direct_path(net, j, i, k)
        ) - sum(
            x[i, j, k, id] for
            j = 1:conf.stations, k = 1:conf.crowd_shippers if csr_direct_path(net, i, j, k)
        ) == 0
    )


    @constraint(
        model,
        select_cs[k = 1:conf.crowd_shippers],
        sum(
            x[i, j, k, p] for
            i = 1:conf.stations, j = 1:conf.stations, p = 1:conf.packages if csr_direct_path(net, i, j, k)
        ) <= 1000z[k]
    )

    @constraint(
        model,
        one_package[
            i = 1:conf.stations,
            j = 1:conf.stations,
            k = 1:conf.crowd_shippers;
            csr_direct_path(net, i, j, k),
        ],
        sum(x[i, j, k, p] for p = 1:conf.packages) <= 1
    )


    @expression(
        model,
        not_before_expr[k = 1:conf.crowd_shippers, p = 1:conf.packages, i = 1:conf.stations],
        sum(arrival[i, k] * x[j, i, k, p] for j ∈ 1:conf.stations if csr_direct_path(net, j, i, k))
    )

    @constraint(
        model,
        not_before[
            k = 1:conf.crowd_shippers,
            p = 1:conf.packages,
            i = 1:conf.stations;
            not_before_expr[k, p, i] ≠ 0.0,
        ],
        t[i, p] >= not_before_expr[k, p, i]
    )


    @expression(
        model,
        not_after_expr[p = 1:conf.packages, i = 1:conf.stations],
        sum(
            arrival[i, k] * x[j, i, k, p] for
            j ∈ 1:conf.stations, k ∈ 1:conf.crowd_shippers if csr_direct_path(net, j, i, k)
        )
    )

    @constraint(
        model,
        not_after[p = 1:conf.packages, i = 1:conf.stations; not_after_expr[p, i] ≠ 0.0],
        t[i, p] <= net.wait[i] + not_after_expr[p, i]
    )

    @expression(
        model,
        positive_time_expr[
            i = 1:conf.stations,
            j = 1:conf.stations,
            p = 1:conf.packages;
            net.adj[i, j] == 1,
        ],
        sum(x[i, j, k, p] for k in 1:conf.crowd_shippers if csr_direct_path(net, i, j, k))
    )

    @constraint(
        model,
        positive_time[
            i = 1:conf.stations,
            j = 1:conf.stations,
            p = 1:conf.packages;
            net.adj[i, j] == 1 && positive_time_expr[i, j, p] ≠ 0,
        ],
        t[i, p] + (net.wait[i] + net.delay[i, j]) * positive_time_expr[i, j, p] <=
        t[j, p] + 1000(1 - positive_time_expr[i, j, p])
    )

    @objective(model, Min, sum(z))

end


function csr_direct_path(net, start, dest, csr)
    net.ct[csr][start, dest] == 1
end
