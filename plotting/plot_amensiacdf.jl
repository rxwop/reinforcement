
using DelimitedFiles, Plots

fascent = readdlm("data\\CDF 1000g fascent.csv", ',', Float64)
ascent = readdlm("data\\CDF 1000g ascent.csv", ',', Float64)
ascent_amnesia = readdlm("data\\CDF 1000g ascent amnesia.csv", ',', Float64)

plot(
    (0:1000) * 70 .+ 1,
    [fascent ascent ascent_amnesia],
    labels = ["fascent" "ascent" "ascent (amnesia)"],
    color = ["0x77c2d8" "0x476eea" "0xB86BE0"],
    xlabel = "Keys Searched",
    ylabel = "Fraction of Runs Solved",
    linewidth = 2.5
)