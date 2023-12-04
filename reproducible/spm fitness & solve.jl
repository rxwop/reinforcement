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
G_set = 500
run_number = 75
N = 50




print("Beginning runs...")

fitness_log = zeros(Float64, G_set + 1)
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
        10, # spawns
        10.; # learning rate (M)
        lineage_habit = :fascent,
        fscore = true
    )

    # initialise solve checker
    solved = false
    solve_num = G_set + 1

    # Add to compute arithmetic mean at the end
    fitness_log[1] += solver.fitness

    for i in 1:G_set
        nextgen!(solver)
        
        fitness_log[i + 1] += solver.fitness
        
        if !solved
            if solver.parent.permutation == inv_target
                solve_num = i
                solved = true
            end
        end
    end

    # Standardise, so f_target = 0.
    fitness_log .-= f_target

    solve_log[run] = solve_num

end

fitness_log /= run_number

solve_percentage = count(!=(G_set + 1), solve_log) / run_number * 100

filter!(!=(G_set + 1), solve_log)

############# EXPORT DATA ###############



println(round(solve_percentage; digits = 1), "% of runs ended in solution")


if !isempty(solve_log)
    println("On average, solved in: ", mean(solve_log), " generations")
end


using Plots
plot(
    fitness_log,
    xlabel = "Generation",
    ylabel = "Mean Relative Fitness"
)
