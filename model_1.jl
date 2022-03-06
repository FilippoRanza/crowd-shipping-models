# Minimize shipment cost
using SparseArrays

using Gurobi
using JuMP

using Configurations
using YAML

@option struct Config
    instance_file::String
    instance_entry::String
end

config = begin
    data = YAML.load_file("config.yml"; dicttype=Dict{String, Any})
    from_dict(Config, data)
end


include("utils.jl")
include("build_instance.jl")
include("load_instance.jl")

inst, traffic, arrival =
    load_instance(config.instance_file, config.instance_entry)


adj = make_adjacent_matrix(inst)

model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "threads", 6)

time_limit(h, m, s) = set_optimizer_attribute(model, "timelimit", h * 3600 + m * 60 + s)

time_limit(0, 5, 0)

conf =
    ModelConfig(size(inst.distance, 1), size(inst.crowd_shippers, 1), size(inst.parcels, 1))
net = NetworkInfo(adj, traffic, inst.wait, inst.distance)

build_model(model, conf, net, inst.parcels, arrival)

optimize!(model)
