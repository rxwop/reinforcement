
using DelimitedFiles, Plots

s5 = readdlm("data\\CDF 1000g 5s ascent.csv", ',', Float64)
s8 = readdlm("data\\CDF 625g 8s ascent.csv", ',', Float64)
s20 = readdlm("data\\CDF 250g 20s ascent.csv", ',', Float64)
s50 = readdlm("data\\CDF 100g 50s ascent.csv", ',', Float64)
s100 = readdlm("data\\CDF 50g 100s ascent.csv", ',', Float64)

color = palette(:viridis, 6)

plot(
    xlabel = "Keys Searched",
    ylabel = "Fraction of Runs Solved"
)

plot!(
    (0:1000) * 5 .+ 1,
    s5,
    label = "5 spawns",
    color = color[1]
)

plot!(
    (0:625) * 8 .+ 1,
    s8,
    label = "8 spawns",
    color = color[2]
)

plot!(
    (0:250) * 20 .+ 1,
    s20,
    label = "20 spawns",
    color = color[3]
)

plot!(
    (0:100) * 50 .+ 1,
    s50,
    label = "50 spawns",
    color = color[4]
)

plot!(
    (0:50) * 100 .+ 1,
    s100,
    label = "100 spawns",
    color = color[5]
)