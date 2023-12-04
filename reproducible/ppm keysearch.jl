cd = @__DIR__

include(joinpath(cd , "..", "ppm reinforcement.jl"))
import Random.randperm, StatsBase.mean

# Load all texts, cut to 2000
using .TxtSamples
tokenise!(orwell)
tokenise!(adfgvx)
tokenise!(jmacro)
texts = [orwell, adfgvx, jmacro]


function avg_keys_searched(G_set::Int, run_number::Int, spawns::Int, M::Float64; bbin::Bool = true, lineage_habit::Symbol = :fascent, texts = texts)

    keys_log = Vector{Int}(undef, run_number)

    for run in 1:run_number

        # Choose random text, note fitness, encrypt with random Substitution
        randtext = copy(rand(texts))[1:2000]
        rand_target = Substitution(randperm(26))
        apply!(rand_target, randtext)

        invert!(rand_target)

        # Choose starting Substitution
        start = invert!(frequency_matched_Substitution(randtext))


        solver = SubstitutionSolve(
            randtext,
            start,
            spawns,
            M; # learning rate (M)
            lineage_habit = lineage_habit,
            bbin = bbin
        )

        # initialise solve checker
        solved = false
        solve_num = 0

        for i in 1:G_set
            nextgen!(solver)
            
            if !solved
                if solver.parent == rand_target
                    solve_num = i
                    solved = true
                end
            end
        end

        if solved
            keys_log[run] = solve_num * spawns + 1
        else
            keys_log[run] = 0
        end

    end

    solve_percentage = count(!=(0), keys_log) / run_number * 100
    keys_searched = mean(filter!(!=(0), keys_log))

    return (keys_searched, solve_percentage)

end

############# EXPORT DATA ###############
