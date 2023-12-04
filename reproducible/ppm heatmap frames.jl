cd = @__DIR__

include(joinpath(cd , "..", "ppm reinforcement.jl"))
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
sample_gen = [1, 10, 20, 30, 50]


print("Beginning...")


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
    lineage_habit = :fascent,
    bbin = true
)


gr(format=:png)
arr_size = solver.ppm.size
p_arr = []

p2 = heatmap(
    perm_matrix(rand_target.mapping),
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
            solver.ppm.arr,
            clims = (0, 1),
            c = cgrad(:matter, rev = true),
            aspect_ratio = :equal,
            legend = :none,
            xlabel = "x",
            ylabel = "y",
            xaxis = false,
            yaxis = false,
            xticks = false,
            yticks = false,
            xlims = (-0.5, arr_size + 0.5),
            ylims = (-0.5, arr_size + 0.5),
            title = "Gen " * string(i)
        ))
    end

end



############# EXPORT DATA ###############
using Measures
p = plot(p_arr..., p2, margin = 5mm, thickness_scaling = 0.7)