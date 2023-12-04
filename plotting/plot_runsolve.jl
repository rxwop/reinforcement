using Plots, DelimitedFiles, Statistics, Plots.PlotMeasures, StatsPlots
include("mycolors.jl")

cd = @__DIR__

names = ["Hill Climb \n(id)", "Hill Climb\n(freq)", "ppM\n(id)", "ppM\n(id, bbin)", "ppM\n(freq)", "ppM\n(freq, bbin)"]
filestems = ["genetic_id", "genetic_f", "ppm_id", "ppm_idbbin", "ppm_f", "ppm_fbbin"]
solves = Vector{Vector{Int}}()
colors = [cocoa[2], viridian[1], cocoa[3], cocoa[4], viridian[3], viridian[4]]

for stem in filestems
    pathend = "data\\sub_bench_solve_" * stem * "_10_fa.csv"
    solve = readdlm(joinpath(cd, "..", pathend), ',', Int)

    push!(solves, vec(solve))
end

p = plot(
)

solves = hcat(solves...)
#solves = permutedims(solves)

boxplot!(p, permutedims(names), solves, legend = false, linecolor = permutedims(colors), fillalpha = 0., markercolor = permutedims(colors), markerstrokewidth = 0, bar_width = 0.5, tickfontsize = 7, thickness_scaling = 1.1)


ylabel!("Solved Generation")

savefig("sub_bench_all_solves_boxp.png")