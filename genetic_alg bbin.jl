push!(LOAD_PATH, raw"C:\Users\robiw\OneDrive\Documents\0_vscode\cipherengine")
using rxciphers


function next_parent_dF(children::Vector{T}, child_dF::Vector{Float64}, lineage_habit = :ascent; parent::T) where T
    if !(length(children) == length(child_dF))
        error("Each child must have an associated delta fitness")
    end

    iter = zip(children, child_dF)

    # Set new parent and fitness
    if lineage_habit == :ascent
        best = first(iter)

        for item in iter # choose best child
            if item[2] > best[2]
                best = item
            end
        end

    elseif lineage_habit == :fascent
        best = (parent, 0.)

        for item in iter # choose best child
            if item[2] > best[2]
                best = item
            end
        end

    elseif lineage_habit == :random
        best = first(iter)

    elseif lineage_habit == :descent
        best = first(iter)

        for item in iter # choose worst child
            if item[2] < best[2]
                best = item
            end
        end

    elseif lineage_habit == :stationary
        best = (parent, 0.)

    else
        error("Invalid lineage habit kwarg")
    end

    return best
end

using StatsBase
function mutate_sample(s::Substitution, indices::CartesianIndices{2, Tuple{Base.OneTo{Int64}, Base.OneTo{Int64}}}, weights::Weights)
    (x1, y2) = sample(indices, weights)
    return (x1, y2)
end

# fitness_log starts from g = 0
# Solves for the inverse substitution
function bench_substitution_solve(start_parent::Substitution, inv_target::Substitution, txt::Txt, gen::Int, spawns::Int; lineage_habit = :ascent; ref_frequencies::Vector{Float64} = monogram_freq)
    parent = start_parent
    parent_f = quadgramlog(parent(txt))

    fitness_log = Vector{Float64}(undef, gen + 1)
    fitness_log[1] = parent_f

    solved = false
    solve_num = gen

    for g in 1:gen
        children = [mutate(parent) for _ in 1:spawns]
        child_fitness = [quadgramlog(child(txt)) for child in children] .- parent_f

        parent, parent_dF = next_parent_dF(children, child_fitness, lineage_habit; parent = parent)
        parent_f += parent_dF

        fitness_log[g + 1] = parent_f

        if !solved
            if parent == inv_target
                solve_num = g
                solved = true
            end
        end

    end

    return fitness_log, solve_num
end
