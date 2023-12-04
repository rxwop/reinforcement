cd = @__DIR__

include(joinpath(cd , "..", "spm reinforcement.jl"))
import Random.randperm
using Plots

# Load all texts, cut to 2000
using .TxtSamples
tokenise!(veryorwell)
# tokenise!(orwell)
# tokenise!(adfgvx)
# tokenise!(jmacro)
# texts = [orwell, adfgvx, jmacro]

# Parameters of data collection
G_set = 100
N = 20
L = 2000
plot_lineage_certainty = false


print("Beginning...")


# Choose random text, note fitness, encrypt with random Permutation
rnd = rand(1:28200 - L)
randtext = copy(veryorwell)[1 + rnd : L + rnd]
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
    fscore = true,
    amnesia_coeff = 0.8,
    amnesia_threshold = 0.95
)


gr(format=:png)
arr_size = solver.spm.size

p1 = heatmap(
    permutedims(solver.spm.arr),
    clims = (0, 1),
    c = cgrad(:deep, rev = true),
    aspect_ratio = :equal,
    xlabel = "x",
    ylabel = "y",
    xaxis = false,
    yaxis = false,
    xticks = false,
    yticks = false,
    xlims = (-0.5, arr_size + 0.5),
    ylims = (-0.5, arr_size + 0.5),
    title = "spM"
)

p2 = heatmap(
    permutedims(perm_matrix(splitsuccessors(inv_target))),
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
    perm_matrix(inv_target, solver.parent.permutation),
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
    p5 = plot([0], [lineage_certainty(solver.spm)],
        xlims = [0, G_set], ylims = [-0.5, 1.5],
        label = "Lineage Certainty" 
    )
else
    p5 = plot([0], [spm_kldiv(solver.spm, splitsuccessors(inv_target))],
    xlims = [0, G_set], ylims = [-0.5, 10.],
    label = "KL Div"
    )
end



using Measures
l = @layout [[p1{0.7w} [p2; p3]]
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
        if solver.parent.permutation == inv_target
            solved.b = true
            scatter!(p[4], [i], [solver.fitness - f_target],
            markershape = :square, linecolor = :black, markercolor = :yellow,
            label = "Solved")
        end
    end

    if solver.fitness > f_target
        println("Spurious key")
    end

    # Plot ppM
    p[1][1][:z] = solver.spm.arr
    p[3][1][:z] = perm_matrix(inv_target, solver.parent.permutation)

    push!(p[4][1][:y], solver.fitness - f_target)
    push!(p[4][1][:x], i)
    if plot_lineage_certainty
        push!(p[5][1][:y], lineage_certainty(solver.spm))
        push!(p[5][1][:x], i)
    else
        push!(p[5][1][:y], spm_kldiv(solver.spm, splitsuccessors(inv_target)))
        push!(p[5][1][:x], i)
    end


end every 1



############# EXPORT DATA ###############

gif(anim, "anim.gif")