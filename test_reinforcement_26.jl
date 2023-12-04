include("ppm reinforcement.jl")

include(raw"C:\Users\robiw\OneDrive\Documents\0_vscode\RxCiphers.jl\test\samples.jl")
using .TxtSamples

W = Alphabet

txt = orwell
tokenise!(txt, W)




S = Affine(11, 8, W)
#S = Substitution(Alphabet)
apply!(S, txt)

@show S




println("Beginning test...")


# using BenchmarkTools
# @btime (PMatrix, cracked) = substitution_solve(txt, 40)
# 30 Hz



# (fitnesses, prefitnesses, solved) = bench_substitution_solve(S, txt, 200, 10, 1.0; lineage_habit = :fascent, bbin = true)

# plot(fitnesses, label = "S fitness")
# plot!(prefitnesses)