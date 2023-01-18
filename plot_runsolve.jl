using Plots, DelimitedFiles, Statistics, Plots.PlotMeasures, StatsPlots
include("mycolors.jl")

names = ["Genetic \n(identity)", "Genetic\n(frequency matched)", "ppM\n(identity, no bbin)", "ppM\n(identity, bbin)", "ppM\n(frequency matched,\nno bbin)", "ppM\n(frequency matched,\nbbin)"]
filestems = ["genetic_id", "genetic_f", "ppm_id", "ppm_idbbin", "ppm_f", "ppm_fbbin"]
solves = Vector{Vector{Int}}()
colors = [cocoa[2], viridian[1], cocoa[3], cocoa[4], viridian[3], viridian[4]]

for stem in filestems
    solve = readdlm("sub_bench_solve_" * stem * "_10_fa.csv", ',', Int)

    push!(solves, vec(solve))
end

p = plot(
    background_color = space,
    foreground_color = cocoa[5],
    margins = 20px
)

solves = hcat(solves...)
#solves = permutedims(solves)

boxplot!(p, permutedims(names), solves, legend = false, linecolor = permutedims(colors), fillcolor = space, markercolor = permutedims(colors), markerstrokewidth = 0, bar_width = 0.5, tickfontsize = 5)


ylabel!("Solved Generation")

savefig("sub_bench_all_solves_boxp.png")