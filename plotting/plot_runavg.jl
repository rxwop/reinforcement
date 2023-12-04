using Plots, DelimitedFiles, Statistics, Plots.PlotMeasures
include("mycolors.jl")

cd = @__DIR__

names = ["Hill Climb (identity)", "Hill Climb (frequency matched)", "ppM (identity, no bbin)", "ppM (identity, bbin)", "ppM (frequency matched, no bbin)", "ppM (frequency matched, bbin)"]
filestems = ["genetic_id", "genetic_f", "ppm_id", "ppm_idbbin", "ppm_f", "ppm_fbbin"]
avgs = Vector{Matrix{Float64}}()
colors = [cocoa[2], viridian[1], cocoa[3], cocoa[4], viridian[3], viridian[4]]

for stem in filestems
    pathend = "data\\sub_bench_fit_" * stem * "_10_fa.csv"
    fitnesses = readdlm(joinpath(cd, "..", pathend), ',', Float64)
    avg = mean(fitnesses; dims = 2)

    push!(avgs, avg)
end

p = plot(
    ylims = (-4, 1),
    xlims = (0, 250)
)

for (n, a, c) in zip(names, avgs, colors)
    a = vec(a)
    plot!(p, a, linewidth = 2.5, label = n, color = c)
end



xlabel!("Generation")
ylabel!("Relative fitness")

savefig("sub_bench_all_avg.png")