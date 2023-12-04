using Plots, DelimitedFiles, Statistics, Plots.PlotMeasures
include("mycolors.jl")

names = ["floored ascent", "ascent", "random", "descent"]
filestems = ["fa", "a", "r", "d"]
avgs = Vector{Matrix{Float64}}()
colors = reverse(habits)

for stem in filestems
    fitnesses = readdlm("data\\sub_bench_fit_ppm_fbbin_10_10_" * stem * ".csv", ',', Float64)
    avg = mean(fitnesses; dims = 2)

    push!(avgs, avg)
end

p = plot(xlims = (0, 500), ylims = (-6, 0.2),
    legend = :right
)

for (n, a, c) in zip(names, avgs, colors)
    plot!(p, a, linewidth = 2, label = n, color = c, subplot = 1)
end

rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])

plot!(p, rectangle(40, 0.5, 20, -0.5), subplot = 1, fillalpha = 0.07, fillcolor = :black, label = false, linealpha = 0.2)

plot!(
    [avgs[1] avgs[2]],
    color = [colors[1] colors[2]],
    xlims = (20, 60),
    ylims = (-0.5, 0.),
    label = false,
    inset = bbox(0.07, -0.02, 150px, 100px, :center),
    subplot = 2
)



ylabel!(p[1], "Relative Fitness")
xlabel!(p[1], "Generation")

savefig("ppm_lhabits.png")