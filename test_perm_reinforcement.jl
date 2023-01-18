include("spm reinforcement.jl")

using JLD2

# @load "jld2/samples.jld2" orwell
# tokenise!(orwell)
@load "jld2/very_orwell.jld2"
tokenise!(very_orwell)

txt = very_orwell[1:2387]


P = Permutation(affine(19, 2, 50), true)
apply!(P, txt)



println("Beginning test...")
# (spM, cracked) = permutation_solve(130, T, 300, 600, 1.0; lineage_habit = :fascent)

(PMatrix, cracked, fitnesses, divergences) = permutation_show_off(P, txt, 200, 100, 1.0; lineage_habit = :fascent, frame_skip = 1)

plot(fitnesses, label = "P fitness")
plot!(divergences / maximum(divergences), label = "spM divergence")

# can you use a higher reinforce_rate?