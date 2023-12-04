cd = @__DIR__

include(joinpath(cd , "..", "ppm reinforcement.jl"))
include(joinpath(cd , "..", "plotting\\mycolors.jl"))
import Random.randperm

# Load all texts, cut to 2000
using .TxtSamples
tokenise!(orwell)
tokenise!(adfgvx)
tokenise!(jmacro)
texts = [orwell, adfgvx, jmacro]

# Parameters of data collection
G_set = 300
run_number = 500




print("Beginning runs...")

fitness_log = Vector{Float64}()
certainty_log = Vector{Float64}()

for run in 1:run_number

    # Choose random text, note fitness, encrypt with random Substitution
    randtext = copy(rand(texts))[1:2000]
    f_target = quadgramlog(randtext)
    rand_target = Substitution(randperm(26))
    apply!(rand_target, randtext)

    invert!(rand_target)

    # Choose starting Substitution
    start = invert!(frequency_matched_Substitution(randtext))


    solver = SubstitutionSolve(
        randtext,
        start,
        10, # spawns
        25.; # learning rate (M)
        lineage_habit = :ascent,
        bbin = false
    )

    for i in 1:G_set
        nextgen!(solver)
        push!(fitness_log, solver.fitness - f_target)
        push!(certainty_log, lineage_certainty(solver.ppm))
        
        if solver.parent == rand_target
            break
        end
    end

end


############# EXPORT DATA ###############

import Statistics.cor, StatsBase.corspearman

println("PMCC: ", cor(fitness_log, certainty_log))
println("Spearman's Rank: ", corspearman(fitness_log, certainty_log))

using Plots

plot_n = 500
using StatsBase
I = sample(1:length(fitness_log), plot_n, replace = false)
scatter(fitness_log[I], certainty_log[I],
    ylabel = "Lineage Certainty",
    xlabel = "Relative Fitness",
    label = "ascent",
    color = habits[3],
    ylims = (-0.04, 1.04)
)