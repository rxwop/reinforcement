cd = @__DIR__

include(joinpath(cd , "..", "spm reinforcement.jl"))
import Random.randperm, Statistics.mean

# Load all texts, cut to 2000
using .TxtSamples
tokenise!(orwell)
tokenise!(adfgvx)
tokenise!(jmacro)
texts = [orwell, adfgvx, jmacro]

# Parameters of data collection
G_set = 1000
run_number = 200
N = 50




print("Beginning runs...")

solve_log = Vector{Int}(undef, run_number)

for run in 1:run_number

    # Choose random text, note fitness, encrypt with random Substitution
    randtext = copy(rand(texts))[1:2000]
    f_target = quadgramlog(randtext)
    rand_target = Permutation(randperm(N))
    apply!(rand_target, randtext)

    inv_target = invert(rand_target).permutation


    solver = PermutationSolve(
        randtext,
        N,
        5, # spawns
        10.; # learning rate (M)
        lineage_habit = :ascent,
        fscore = true,
    )

    # initialise solve checker
    solved = false
    solve_num = G_set + 1

    for i in 1:G_set
        nextgen!(solver)
    
        
        if !solved
            if solver.parent.permutation == inv_target
                solve_num = i
                solved = true
            end
        end

        if solved
            break
        end
    end

    solve_log[run] = solve_num

end

############# EXPORT DATA ###############

using StatsBase
c = ecdf(solve_log)

Y = c.(0:G_set)

using DelimitedFiles
writedlm("data\\CDF 1000g 50s ascent.csv", Y)

# using Plots
# plot(
#     0:G_set, x -> c(x),
# )



## CDF 500g - N = 50, L = 2000, 10 spawns
## CDF 1000g - N = 100, L = 2000, 70 spawns
