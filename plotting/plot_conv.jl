using Plots, DelimitedFiles
divergences = readdlm(raw"data\sub_bench_div_ppm_fbbin_10_25_fa_100.csv", ',', Float64)
variances = readdlm(raw"data\sub_bench_var_ppm_fbbin_10_25_fa_100.csv", ',', Float64)

include("mycolors.jl")

using Plots.PlotMeasures

p = plot(thickness_scaling = 1.15
)



using Statistics
avgdiv = mean(divergences; dims = 2)
avgvar = mean(variances; dims = 2)

plot!(p, 0:100, avgdiv, linewidth = 2, label = "Mean KL Divergence")

plot!(p, 0:100, avgvar, linewidth = 2, label = "Mean Col. Variance Sum")

xlabel!("Generation")


savefig("convergence ppm_fbbin_10_25_fa_100.png")