using Plots, DelimitedFiles, Statistics, Plots.PlotMeasures
include("mycolors.jl")

names = ["floored ascent", "ascent", "random", "descent"]
filestems = ["fa", "a", "r", "d"]
avgs = Vector{Matrix{Float64}}()
colors = reverse!(cocoa[2:5])

for stem in filestems
    fitnesses = readdlm("sub_bench_fit_genetic_f_10_" * stem * ".csv", ',', Float64)
    avg = mean(fitnesses; dims = 2)

    push!(avgs, avg)
end

p = plot(xlims = (0, 500), ylims = (-6, 0.2),
    background_color = space,
    foreground_color = cocoa[5],
    margins = 20px,
    legend = :right
)

for (n, a, c) in zip(names, avgs, colors)
    plot!(p, a, linewidth = 2, label = n, color = c)
end



ylabel!("Relative Fitness")
xlabel!("Generation")

savefig("genetic_lhabits.png")