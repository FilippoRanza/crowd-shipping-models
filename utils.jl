
function print_result(model, pkgs, stats, csr, adj, trans)
    if termination_status(model) == OPTIMAL
        for p in 1:pkgs
            println("Package $p")
            for k in 1:csr
                output = false
                for i in 1:stats, j in 1:stats 
                    if adj[i, j] == 1 && trans[k][i, j] == 1
                        if value(model[:x][i, j, k, p]) == 1.
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

function print_element(model::Model, key::Symbol)
    element = model[key]
    for e in element
        println(e)
    end
end