using Plots, DelimitedFiles
fitnesses = readdlm("sub_bench_fit_ppm_fbbin_10_10_d_5000.csv", ',', Float64)

include("mycolors.jl")

using Plots.PlotMeasures

p = plot(xlims = (0, 5000), ylims = (-Inf, 1),
    background_color = space,
    foreground_color = cocoa[5],
    margins = 20px
)

min = minimum(fitnesses; dims = 2)
max = maximum(fitnesses; dims = 2)


using Statistics
avg = mean(fitnesses; dims = 2)

plot!(p, 0:5000, min, linewidth = 0, label = false, line_z = min, color = cgrad(viridian), colorbar = false)
plot!(p, 0:5000, max, fillrange = min, fillalpha = 0.35, linewidth = 0, label = false, fill_z = avg, color = cgrad(viridian), colorbar = false)

plot!(p, 0:5000, avg, linewidth = 3, label = "Average (mean)", line_z = avg, color = cgrad(viridian), colorbar = false)



xlabel!("Generation")
ylabel!("Relative fitness")

savefig("ppm_descent_5000.png")