using RxCiphers, Plots

###### SETUP CIPHERTEXT

include("ppm reinforcement.jl")

target = Atbash(676)
text = tokenise!(TxtSamples.veryorwell)[1:10000]
nchar!(text, 2)

ctext = target(text)


###### 

G_set = 1000
plot_lineage_certainty = true


####### 

print("Beginning...")

f_target = quadgramlog(nchar(text, 1))

# Choose starting Substitution
start = Substitution(676)


solver = BiSubSolve(
    ctext,
    start,
    200, # spawns
    10.; # learning rate (M)
    lineage_habit = :ascent,
    bbin = false
)


###### PLOTTING #####################

invert!(target)


arr_size = solver.ppm.size

gr(format=:png)
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
    perm_matrix(target.mapping),
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
    perm_matrix(target.mapping, solver.parent.mapping),
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
    p5 = plot([0], [ppm_kldiv(solver.ppm, target.mapping)],
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
    
    println("Generation ", i, " beginning...")
    nextgen!(solver)

    if !solved.b
        if solver.parent == target
            solved.b = true
            scatter!(p[4], [i], [solver.fitness - f_target],
            markershape = :square, linecolor = :black, markercolor = :yellow,
            label = "Solved")
        end
    end

    # Plot ppM
    p[1][1][:z] = solver.ppm.arr
    p[3][1][:z] = perm_matrix(target.mapping, solver.parent.mapping)


    push!(p[4][1][:y], solver.fitness - f_target)
    push!(p[4][1][:x], i)
    if plot_lineage_certainty
        push!(p[5][1][:y], lineage_certainty(solver.ppm))
        push!(p[5][1][:x], i)
    else
        push!(p[5][1][:y], ppm_kldiv(solver.ppm, target.mapping))
        push!(p[5][1][:x], i)
    end


end every 10



############# EXPORT DATA ###############

gif(anim, "anim.gif")