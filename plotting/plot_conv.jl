using Plots, DelimitedFiles
divergences = readdlm(raw"data\sub_bench_div_ppm_fbbin_10_25_fa_100.csv", ',', Float64)
variances = readdlm(raw"data\sub_bench_var_ppm_fbbin_10_25_fa_100.csv", ',', Float64)

include("mycolors.jl")

using Plots.PlotMeasures

p = plot(xlims = (0, 100), ylims = (0, Inf),
background_color = space,
foreground_color = cocoa[5],
margins = 20px
)



using Statistics
avgdiv = mean(divergences; dims = 2)
avgvar = mean(variances; dims = 2)

plot!(p, 0:100, avgdiv, linewidth = 3, label = "Mean KL-Divergence", line_z = avgdiv, color = cgrad(viridian[1:2]), colorbar = false)

plot!(p, 0:100, avgvar, linewidth = 3, label = "Mean variance", line_z = avgvar, color = cgrad(viridian[3:4]), colorbar = false)

xlabel!("Generation")


savefig("convergence ppm_fbbin_10_25_fa_100.png")