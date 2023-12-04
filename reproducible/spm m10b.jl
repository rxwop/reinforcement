cd = @__DIR__

using RxCiphers
text = Txt(read(raw"reproducible\\m10b.txt", String))
tokenise!(text)

include(joinpath(cd , "..", "spm reinforcement.jl"))
using Plots

# Parameters of data collection
G_set = 600
plot_lineage_certainty = true


m = collect(1:130)
m = reshape(m, (5, :))
for (i, j) in enumerate([22, 18, 15, 14, 8])
    m[i, :] = circshift(m[i, :], - j)
end
m = m[[4,2,3,1,5], :]
target = invPermutation(vec(m))

plain = target(text)
f_target = quadgramlog(plain)

print("Beginning...")



solver = PermutationSolve(
    text,
    130,
    150, # spawns
    15.; # learning rate (M)
    lineage_habit = :ascent,
    fscore = true,
    amnesia_coeff = 0.8,
    amnesia_threshold = 0.95
)


gr(format=:png)
arr_size = solver.spm.size


p1 = heatmap(
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
    xlims = (-0.5, arr_size + 0.5),
    ylims = (-0.5, arr_size + 0.5),
    title = "spM"
)

p2 = heatmap(
    permutedims(perm_matrix(splitsuccessors(target.permutation))),
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
    perm_matrix(target.permutation, solver.parent.permutation),
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
    p5 = plot([0], [spm_kldiv(solver.spm, splitsuccessors(target.permutation))],
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
        if solver.parent.permutation == target.permutation
            solved.b = true
            scatter!(p[4], [i], [solver.fitness - f_target],
            markershape = :square, linecolor = :black, markercolor = :yellow,
            label = "Solved")
        end
    end

    # Plot ppM
    p[1][1][:z] = solver.spm.arr
    p[3][1][:z] = perm_matrix(target.permutation, solver.parent.permutation)

    push!(p[4][1][:y], solver.fitness - f_target)
    push!(p[4][1][:x], i)
    if plot_lineage_certainty
        push!(p[5][1][:y], lineage_certainty(solver.spm))
        push!(p[5][1][:x], i)
    else
        push!(p[5][1][:y], spm_kldiv(solver.spm, splitsuccessors(target.permutation)))
        push!(p[5][1][:x], i)
    end


end every 20




############# EXPORT DATA ###############

gif(anim, "anim.gif")