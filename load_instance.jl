


using HDF5


struct LoadInstace
    parcels::Any
    crowd_shippers::Any
    distance::Any
    wait::Any
end

function load_instance(file_name, entry_name)
    h5open(file_name) do file
        load(file[entry_name])
    end
end

function make_adjacent_matrix(li::LoadInstace)
    sz = size(li.distance, 1)
    adj = ones(sz, sz)
    for i = 1:sz
        adj[i, i] = 0
    end
    adj
end

function load(data)
    parcels = read(data["parcels"])
    csrs = read(data["csrs"])
    distance = read(data["delay"])
    wait = read(data["wait"])
    li = LoadInstace(parcels, csrs, distance, wait)
    T = make_sparse_csr_array(li)
    A = make_arrive_matrix(li)
    li, T, A
end


function make_sparse_csr_array(li::LoadInstace)
    stations = size(li.distance, 1)
    csr = size(li.crowd_shippers, 1)
    T = [spzeros(Int8, stations, stations) for i = 1:csr]
    for i = 1:csr
        s, d = li.crowd_shippers[i, 2:3]
        T[i][s, d] = 1
    end
    T
end


function make_arrive_matrix(li::LoadInstace)
    c = size(li.crowd_shippers, 1)
    s = size(li.distance, 1)
    A = Inf * ones(s, c)
    for (i, (w, s, d)) in enumerate(eachrow(li.crowd_shippers))
        A[s, i] = w
        A[d, i] = w + li.wait[s] + li.distance[s, d]
    end
    A
end
