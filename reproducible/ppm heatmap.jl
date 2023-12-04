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


print("Beginning...")


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


gr(format=:png)
arr_size = solver.ppm.size

p1 = heatmap(
    permutedims(solver.ppm.arr),
    clims = (0, 1),
    c = cgrad(:matter, rev = true),
    aspect_ratio = :equal,
    xlabel = "x",
    ylabel = "y",
    xaxis = false,
    yaxis = false,
    xticks = false,
    yticks = false,
    xlims = (-0.5, arr_size + 0.5),
    ylims = (-0.5, arr_size + 0.5),
    title = "ppM"
)

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

using Measures
l = @layout [p1{0.7w} p2]
p = plot(p1, p2, layout = l, margin = 3mm)

anim = @animate for i in 1:G_set
    nextgen!(solver)

    # Plot ppM
    p[1][1][:z] = solver.ppm.arr

end every 1



############# EXPORT DATA ###############

gif(anim, "anim.gif")