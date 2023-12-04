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
sample_gen = [1, 5, 10, 15, 20, 50]


print("Beginning...")


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
    lineage_habit = :fascent,
    fscore = true
)


gr(format=:png)
arr_size = solver.spm.size
p_arr = []



for i in 1:G_set
    nextgen!(solver)

    if i in sample_gen
        # Plot ppM
        println(inv_target, solver.parent.permutation)

        push!(p_arr, heatmap(
            perm_matrix(inv_target, solver.parent.permutation),
            clims = (0, 1),
            c = cgrad(["0x131314", "0xFAFAFA"]),
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
p = plot(p_arr..., margin = 5mm, thickness_scaling = 0.7)