#include("genetic_alg.jl")
include("ppm reinforcement.jl")


gen_set = 75
run_number = 500



include(raw"C:\Users\robiw\OneDrive\Documents\0_vscode\RxCiphers.jl\test\samples.jl")
using .TxtSamples

tokenise!(orwell)
tokenise!(adfgvx)
tokenise!(jmacro)

txts = [orwell, adfgvx, jmacro]


print("Beginning benchmark...")

fitnesses = Vector{Float64}(undef, run_number)
solve_lengths = Vector{Int}(undef, run_number)
divergences = Vector{Float64}(undef, run_number)


using Random
for i in 1:run_number
    randtxt = copy(rand(txts))[1:2000]
    f_target = quadgramlog(randtxt)
    rand_target = Substitution(randperm(26))
    apply!(rand_target, randtxt)

    start = invert!(frequency_matched_Substitution(randtxt))


    ## RUN
    #fit, gen = bench_substitution_solve(start, invert!(rand_target), randtxt, gen_set, 10; lineage_habit = :descent)
    fit, div, gen = enddata_substitution_solve(start, rand_target, invert(rand_target), randtxt, gen_set, 10, Inf;
    lineage_habit = :fascent, bbin = true)


    fitnesses[i] = fit - f_target
    solve_lengths[i] = gen
    divergences[i] = div

    println(i)
end


using Statistics
solves_only = filter(!=(gen_set + 1), solve_lengths)
avg_solve = mean(solves_only)
stdev_solve = stdm(solves_only, avg_solve)

frac_solve = count(!=(gen_set + 1), solve_lengths) / run_number * 100


println("Solved in average: $(avg_solve) generations (σ = $(stdev_solve))")
println("$frac_solve % solved")