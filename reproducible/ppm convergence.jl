cd = @__DIR__

include(joinpath(cd , "..", "ppm reinforcement.jl"))
import Random.randperm, Statistics.var

# Load all texts, cut to 2000
using .TxtSamples
tokenise!(orwell)
tokenise!(adfgvx)
tokenise!(jmacro)
texts = [orwell, adfgvx, jmacro]

# Parameters of data collection
G_set = 100
run_number = 150




print("Beginning runs...")

div_log = zeros(Float64, G_set + 1)
var_log = zeros(Float64, G_set + 1)

for run in 1:run_number

    # Choose random text, note fitness, encrypt with random Substitution
    randtext = copy(rand(texts))[1:2000]
    f_target = quadgramlog(randtext)
    rand_target = Substitution(randperm(26))
    apply!(rand_target, randtext)

    # Choose starting Substitution
    start = invert!(frequency_matched_Substitution(randtext))


    solver = SubstitutionSolve(
        randtext,
        start,
        10, # spawns
        25.; # learning rate (M)
        lineage_habit = :fascent,
        bbin = true
    )

    # Add to compute arithmetic mean at the end
    div_log[1] += ppm_kldiv(solver.ppm, rand_target.mapping)
    var_log[1] += sum(var(solver.ppm.arr; dims = 2))

    for i in 1:G_set
        nextgen!(solver)
        div_log[i + 1] += ppm_kldiv(solver.ppm, rand_target.mapping)
        var_log[i + 1] += sum(var(solver.ppm.arr; dims = 2))
    end


end

div_log /= run_number
var_log /= run_number

############# EXPORT DATA ###############

using Plots
plot(
    [div_log var_log],
    label = ["Mean KL Divergence" "Mean Col. Variance Sum"]
)