# Minimize shipment cost
using SparseArrays

using Gurobi
using JuMP

using Configurations
using YAML

using CSV
using DataFrames

@option struct Config
    instance_file::String
    instance_entry::String
    solution_dir::String
    stats_file::String
    result_file::String
    start_id::Int64
    end_id::Int64
end

config = begin
    data = YAML.load_file("config.yml"; dicttype = Dict{String,Any})
    from_dict(Config, data)
end



include("utils.jl")
include("build_instance.jl")
include("load_instance.jl")

results =
    DataFrame(time = Float64[], nodes = Float64[], has_sol = Bool[], sol_value = Float64[])

function make_dataset_entry(m::Model)
    sol_time = solve_time(m)
    nodes = node_count(m)
    has_sol = has_values(m)
    sol_value = if has_sol
        objective_value(m)
    else
        0
    end
    [sol_time, nodes, has_sol, sol_value]
end

function update_dataset!(df::DataFrame, m::Model)
    entry = make_dataset_entry(m)
    push!(df, entry)
end

function save_result(sol_path, m::Model, id)
    x = value.(m[:x])
    open(joinpath(sol_path, "solution-$id.txt"), "w") do file
        for (k, v) in pairs(x.data)
            if v == 1
                println(file, k)
            end
        end
    end
end

function avg_constr_len(model, constr)
    constr_ref = model[constr]
    number = length(constr_ref)
    total = 0
    for c_ref in constr_ref
        c = constraint_object(c_ref)
        total += length(c.func.terms)
    end
    total / number
end

function add_stat_df_entry!(df::DataFrame, m::Model, csrs, prls)
    entry = [
        csrs,
        prls,
        length(m[:x]),
        length(m[:t]),
        avg_constr_len(m, :start_pkg),
        avg_constr_len(m, :arrive_pkg),
        avg_constr_len(m, :flow),
        avg_constr_len(m, :not_before),
        avg_constr_len(m, :not_after),
        avg_constr_len(m, :positive_time),
    ]
    push!(df, entry)
end

stats_df = DataFrame(
    csr = Int64[],
    prcl = Int64[],
    x_vars = Int64[],
    t_vars = Int64[],
    start_pkg = Float64[],
    arrive_pkg = Float64[],
    flow = Float64[],
    not_before = Float64[],
    not_after = Float64[],
    positive_time = Float64[],
)

model = for i ∈ config.start_id:config.end_id
    entry = "$(config.instance_entry)-$i"
    println(entry)
    inst, traffic, arrival = load_instance(config.instance_file, entry)
    adj = make_adjacent_matrix(inst)

    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "threads", 6)


    time_limit(h, m, s) =
        set_optimizer_attribute(model, "timelimit", h * 3600 + m * 60 + s)

    time_limit(0, 15, 0)

    n_stations = size(inst.distance, 1)
    n_csrs = size(inst.crowd_shippers, 1)
    n_parcels = size(inst.parcels, 1)

    conf = ModelConfig(n_stations, n_csrs, n_parcels)
    net = NetworkInfo(adj, traffic, inst.wait, inst.distance)

    build_model(model, conf, net, inst.parcels, arrival)
    add_stat_df_entry!(stats_df, model, n_csrs, n_parcels)

    optimize!(model)

    update_dataset!(results, model)
    if has_values(model)
        save_result(config.solution_dir, model, i)
    end
end

open(config.stats_file, "w") do file
    CSV.write(file, stats_df)
end


open(config.result_file, "w") do file
    CSV.write(file, results)
end
