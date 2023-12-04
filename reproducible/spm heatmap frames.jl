cd = @__DIR__

include(joinpath(cd , "..", "spm reinforcement.jl"))
import Random.randperm
using Plots

# Load all texts, cut to 2000
using .TxtSamples
tokenise!(orwell)
tokenise!(adfgvx)
tokenise!(jmacro)
texts = [orwell, adfgvx, jmacro]

# Parameters of data collection
G_set = 60
N = 20
sample_gen = [1, 3, 5, 8, 10]


print("Beginning...")


# Choose random text, note fitness, encrypt with random Permutation
randtext = copy(rand(texts))[1:2000]
f_target = quadgramlog(randtext)
rand_target = Permutation(randperm(26))
apply!(rand_target, randtext)

inv_target = invperm(rand_target.permutation)


solver = PermutationSolve(
    randtext,
    N,
    10, # spawns
    25.; # learning rate (M)
    lineage_habit = :fascent,
    fscore = true
)


gr(format=:png)
arr_size = solver.spm.size
p_arr = []

p2 = heatmap(
    perm_matrix(splitsuccessors(inv_target)),
    clims = (0, 1),
    c = cgrad(["0x131314", "0xFAFAFA"]),
    legend = :none,
    aspect_ratio = :equal,
    xlabel = "x",
    ylabel = "y",
    xaxis = false,
    yaxis = false,
    xticks = false,
    yticks = false,
    xlims = (-0.5, arr_size + 0.5),
    ylims = (-0.5, arr_size + 0.5),
    title = "Target"
)



for i in 1:G_set
    nextgen!(solver)

    if i in sample_gen
        # Plot ppM
        push!(p_arr, heatmap(
            solver.spm.arr,
            clims = (0, 1),
            c = cgrad(:deep, rev = true),
            aspect_ratio = :equal,
            xlabel = "x",
            ylabel = "y",
            xaxis = false,
            yaxis = false,
            xticks = false,
            yticks = false,
            colorbar = false,
            xlims = (-0.5, arr_size + 0.5),
            ylims = (-0.5, arr_size + 0.5),
            title = "Gen " * string(i)
        ))
    end

end



############# EXPORT DATA ###############
using Measures
p = plot(p_arr..., p2, margin = 5mm, thickness_scaling = 0.7)