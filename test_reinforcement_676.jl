include("reinforcement2.jl")

using JLD2

@load "jld2/samples.jld2" orwell
@load "jld2/bigram_frequencies_vec.jld2" bigram_freq_vec

W = Alphabet ^ 2

txt = orwell
tokenise!(txt, W)




S = Affine(3, 5, W)
#S = Substitution(Alphabet)
apply!(S, txt)

@show S




println("Beginning test...")


# using BenchmarkTools
# @btime (PMatrix, cracked) = linear_reinforcement(txt, 100, 10, Choice_Weights, quadgramlog, eng, 3.0; lineage_habit = "floored ascent")



(PMatrix, cracked, fitnesses, divergences) = debug_cw_substitution_solve(invert(S), txt, 800, 50000, bibigramlog_arr, bigram_freq_vec, 0.5; lineage_habit = :ascent, frame_skip = 100)

plot(fitnesses, label = "S fitness")
plot!(divergences / maximum(divergences), label = "ppM divergence")