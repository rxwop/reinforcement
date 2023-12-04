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
plot_lineage_certainty = false


print("Beginning...")


# Choose random text, note fitness, encrypt with random Substitution
randtext = copy(rand(texts))[1:2000]
f_target = quadgramlog(randtext)
rand_target = Substitution(randperm(26))
apply!(rand_target, randtext)

inv_rand_target = invert(rand_target)


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

arr_size = solver.ppm.size

gr(format=:png)
p1 = heatmap(
    permutedims(getprob!(solver.ppm)),
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

p3 = heatmap(
    perm_matrix(inv_rand_target.mapping, solver.parent.mapping),
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
    title = "Permutation to Soln."
)

p4 = plot([0], [solver.fitness - f_target],
    xlims = [0, G_set], ylims = [-6., 0.5],
    label = "Fitness"
)


if plot_lineage_certainty
    p5 = plot([0], [lineage_certainty(solver.ppm)],
        xlims = [0, G_set], ylims = [-0.5, 1.5],
        label = "Lineage Certainty" 
    )
else
    p5 = plot([0], [ppm_kldiv(solver.ppm, rand_target.mapping)],
        xlims = [0, G_set], ylims = [-0.5, 10.],
        label = "KL Div"
    )
end



using Measures
l = @layout [[p1{0.7w} [p2
p3]]
grid(1, 2){0.3h}]
p = plot(p1, p2, p3, p4, p5, thickness_scaling = 0.7, margin = 3mm, layout = l)

# initialise solve checker
mutable struct Solved
    b::Bool
    Solved() = new(false)
end
solved = Solved()

anim = @animate for i in 1:G_set
    
    nextgen!(solver)

    if !solved.b
        if solver.parent == inv_rand_target
            solved.b = true
            scatter!(p[4], [i], [solver.fitness - f_target],
            markershape = :square, linecolor = :black, markercolor = :yellow,
            label = "Solved")
        end
    end

    # Plot ppM
    p[1][1][:z] = getprob!(solver.ppm)
    p[3][1][:z] = perm_matrix(inv_rand_target.mapping, solver.parent.mapping)


    push!(p[4][1][:y], solver.fitness - f_target)
    push!(p[4][1][:x], i)
    if plot_lineage_certainty
        push!(p[5][1][:y], lineage_certainty(solver.ppm))
        push!(p[5][1][:x], i)
    else
        push!(p[5][1][:y], ppm_kldiv(solver.ppm, rand_target.mapping))
        push!(p[5][1][:x], i)
    end


end every 1



############# EXPORT DATA ###############

gif(anim, "anim.gif")