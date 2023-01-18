using Plots, DelimitedFiles
fitnesses = readdlm("sub_bench_fit_ppm_fbbin_10_10_d_5000.csv", ',', Float64)

include("mycolors.jl")

using Plots.PlotMeasures

p = plot(xlims = (0, 5000), ylims = (-Inf, 1),
    background_color = space,
    foreground_color = cocoa[5],
    margins = 20px
)

for col in eachcol(fitnesses)
    plot!(p, 0:5000, col, linealpha = 0.025, linewidth = 3, label = false, line_z = col, color = cgrad(viridian), colorbar = false)
end

using Statistics
avg = mean(fitnesses; dims = 2)

plot!(p, 0:5000, avg, linewidth = 2, label = "Average (mean)", line_z = avg, color = cgrad(viridian), colorbar = false)



xlabel!("Generation")
ylabel!("Relative fitness")

savefig("ppm_descent_5000.png")