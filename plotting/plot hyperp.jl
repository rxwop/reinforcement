using Plots, JLD2

@load "data\\hyperp narrow1_fbbin.jld2"

include("plotting\\mycolors.jl")

using Plots.PlotMeasures

p = plot(
)



using Statistics
solves_only = filter.(!=(76), solve_lengths)
avgsolve = mean.(solves_only)
stdsolve = std.(solves_only)
percentsolve = count.(!=(76), solve_lengths) / 75
avgdiv = mean.(divergences)

iter_rate = 0:2:20
iter_spawns = 1:5:35

searched = avgsolve .* permutedims(collect(iter_spawns)) .+ 1

data = percentsolve

heatmap!(p, iter_spawns, iter_rate, data, color = cgrad(reverse(viridian)), aspect_ratio = :equal, colorbar = false,
xticks = iter_spawns, yticks = iter_rate, grid = false,
xlims = (-1.5, 33.5), ylims = (-1, 21))

for i in 1:11
    for j in 1:7
        annotate!(p,
            iter_spawns[j],
            iter_rate[i],
            text(round(data[i, j]; digits = 2); pointsize = 8),
        )
    end
end

ylabel!("Update Rate")
xlabel!("Spawn Number")


savefig("hyperp narrow 1 fbbin percentsolve.png")