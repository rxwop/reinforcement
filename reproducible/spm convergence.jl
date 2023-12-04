cd = @__DIR__

include(joinpath(cd , "..", "spm reinforcement.jl"))
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
N = 20




print("Beginning runs...")

div_log = zeros(Float64, G_set + 1)
div_num = fill(run_number, G_set + 1)
var_log = zeros(Float64, G_set + 1)

for run in 1:run_number

    # Choose random text, note fitness, encrypt with random Permutation
    randtext = copy(rand(texts))[1:2000]
    f_target = quadgramlog(randtext)
    rand_target = Permutation(randperm(N))
    apply!(rand_target, randtext)

    inv_target = invperm(rand_target.permutation)


    solver = PermutationSolve(
        randtext,
        N,
        10, # spawns
        25.; # learning rate (M)
        lineage_habit = :ascent,
        fscore = true
    )

    # Add to compute arithmetic mean at the end
    d = spm_kldiv(solver.spm, splitsuccessors(inv_target))
    if !isinf(d)
        div_log[1] += d
    else
        div_num[1] -= 1
    end
    var_log[1] += sum(var(solver.spm.arr; dims = 2))

    for i in 1:G_set
        nextgen!(solver)
        d = spm_kldiv(solver.spm, splitsuccessors(inv_target))
        if !isinf(d)
            div_log[i + 1] += d
        else
            div_num[i + 1] -= 1
        end
        var_log[i + 1] += sum(var(solver.spm.arr; dims = 2))
    end


end

div_log ./= div_num
var_log /= run_number

############# EXPORT DATA ###############

using Plots
plot(
    [div_log var_log],
    label = ["Mean KL Divergence" "Mean Col. Variance Sum"],
    xlabel = "Generation",
    thickness_scaling = 1.1
)

plot!(div_num / run_number, label = "Fraction of Useable Runs")