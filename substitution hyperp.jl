
include("ppm reinforcement.jl")


gen_set = 75
run_number = 75



include(raw"C:\Users\robiw\OneDrive\Documents\0_vscode\RxCiphers.jl\test\samples.jl")
using .TxtSamples

tokenise!(orwell)
tokenise!(adfgvx)
tokenise!(jmacro)

txts = [orwell, adfgvx, jmacro]


print("Beginning benchmark...")


using IterTools
iter_rate = 1:10:50
iter_spawns = 1:10:100
iter = IterTools.product(iter_rate, iter_spawns)

m_size = size(iter)

fitnesses = Matrix{Vector{Float64}}(undef, m_size)
solve_lengths = Matrix{Vector{Int}}(undef, m_size)
divergences = Matrix{Vector{Float64}}(undef, m_size)


using Random

for (i, (rt, sp)) in enumerate(iter)
    run_fitnesses = Vector{Float64}(undef, run_number)
    run_solves = Vector{Int}(undef, run_number)
    run_divergences = Vector{Float64}(undef, run_number)

    for n in 1:run_number
        randtxt = copy(rand(txts))[1:2000]
        f_target = quadgramlog(randtxt)
        rand_target = Substitution(randperm(26))
        apply!(rand_target, randtxt)

        start = invert!(frequency_matched_Substitution(randtxt))


        ## RUN
        fit, div, solve = enddata_substitution_solve(start, rand_target, invert(rand_target), randtxt, gen_set, sp, float(rt); lineage_habit = :fascent, bbin = true)


        run_fitnesses[n] = fit - f_target
        run_solves[n] = solve
        run_divergences[n] = div
    end

    fitnesses[i] = run_fitnesses
    solve_lengths[i] = run_solves
    divergences[i] = run_divergences
    println(i)
end



# genetic _ start _ spawns _ lineage_habit
name = "_fbbin"

using JLD2
jldsave("hyperp" * name * ".jld2"; fitnesses, solve_lengths, divergences)