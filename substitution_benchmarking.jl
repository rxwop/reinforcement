#include("genetic_alg.jl")
include("ppm reinforcement.jl")


gen_set = 5000
run_number = 150



using JLD2
@load "jld2/samples.jld2" orwell
@load "jld2/moresamples.jld2"

tokenise!(orwell)
tokenise!(adfgvx)
tokenise!(jmacro)

txts = [orwell, adfgvx, jmacro]


print("Beginning benchmark...")

fitnesses = Vector{Vector{Float64}}(undef, run_number)
solve_lengths = Vector{Int}(undef, run_number)


using Random
for i in 1:run_number
    randtxt = copy(rand(txts))[1:2000]
    f_target = quadgramlog(randtxt)
    rand_target = Substitution(randperm(26))
    apply!(rand_target, randtxt)

    start = invert!(frequency_matched_Substitution(randtxt))


    ## RUN
    #fit, gen = bench_substitution_solve(start, invert!(rand_target), randtxt, gen_set, 10; lineage_habit = :descent)
    fit, gen = bench_substitution_solve(start, invert!(rand_target), randtxt, gen_set, 10, 10.0; lineage_habit = :descent, bbin = true)


    fitnesses[i] = fit .- f_target
    solve_lengths[i] = gen

    println(i)
end


fitnesses = hcat(fitnesses...)

avg_fitness = sum(fitnesses, dims = 2) / run_number


using Statistics
avg_solve = mean(solve_lengths)
stdev_solve = stdm(solve_lengths, avg_solve)

if any(solve_lengths .== gen_set)
    println("Warning: not all runs resulted in solution")
end



# genetic _ start _ spawns _ lineage_habit
name = "ppm_fbbin_10_10_d_5000"

using DelimitedFiles
writedlm("sub_bench_fit_" * name * ".csv", fitnesses, ',')
writedlm("sub_bench_solve_" * name * ".csv", solve_lengths, ",")
println("Solved in average: $(avg_solve) generations (σ = $(stdev_solve))")



using Plots
plot(avg_fitness)