using Plots, DelimitedFiles, Statistics, Plots.PlotMeasures
include("mycolors.jl")

names = ["Genetic (identity)", "Genetic (frequency matched)", "ppM (identity, no bbin)", "ppM (identity, bbin)", "ppM (frequency matched, no bbin)", "ppM (frequency matched, bbin)"]
filestems = ["genetic_id", "genetic_f", "ppm_id", "ppm_idbbin", "ppm_f", "ppm_fbbin"]
avgs = Vector{Matrix{Float64}}()
colors = [cocoa[2], viridian[1], cocoa[3], cocoa[4], viridian[3], viridian[4]]

for stem in filestems
    fitnesses = readdlm("sub_bench_fit_" * stem * "_10_fa.csv", ',', Float64)
    avg = mean(fitnesses; dims = 2)

    push!(avgs, avg)
end

p = plot(xlims = (-2, 0), ylims = (0, 0.2),
    background_color = space,
    foreground_color = cocoa[5],
    margins = 20px,
    legend = :topright
)

for (n, a, c) in zip(names, avgs, colors)
    a = vec(a)
    da = diff(a)
    pop!(a)
    plot!(p, a, da, linewidth = 1, label = n, color = c)
end



ylabel!("Delta Fitness (per generation)")
xlabel!("Relative fitness")

savefig("sub_bench_all_avg_diff_cropped.png")