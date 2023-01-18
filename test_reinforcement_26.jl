include("ppm reinforcement.jl")

using JLD2

@load "jld2/samples.jld2" orwell

W = Alphabet

txt = orwell
tokenise!(txt, W)




S = Affine(11, 8, W)
#S = Substitution(Alphabet)
apply!(S, txt)

@show S




println("Beginning test...")


# using BenchmarkTools
# @btime (PMatrix, cracked) = linear_reinforcement(txt, 100, 10, Choice_Weights, quadgramlog, eng, 3.0; lineage_habit = "floored ascent")



(PMatrix, cracked, fitnesses, divergences) = substitution_show_off(S, txt, 200, 10, 1.0; lineage_habit = :fascent)

plot(fitnesses, label = "S fitness")
plot!(divergences / maximum(divergences), label = "ppM divergence")