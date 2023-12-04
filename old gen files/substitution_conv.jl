
include("ppm reinforcement.jl")


gen_set = 100
run_number = 500



include(raw"C:\Users\robiw\OneDrive\Documents\0_vscode\RxCiphers.jl\test\samples.jl")
using .TxtSamples

tokenise!(orwell)
tokenise!(adfgvx)
tokenise!(jmacro)

txts = [orwell, adfgvx, jmacro]


print("Beginning benchmark...")

divergences = Vector{Vector{Float64}}(undef, run_number)
variances = Vector{Vector{Float64}}(undef, run_number)


using Random
for i in 1:run_number
    randtxt = copy(rand(txts))[1:2000]
    f_target = quadgramlog(randtxt)
    rand_target = Substitution(randperm(26))
    apply!(rand_target, randtxt)

    start = invert!(frequency_matched_Substitution(randtxt))


    ## RUN
    div, var = conv_substitution_solve(start, rand_target, randtxt, gen_set, 10, 25.0; lineage_habit = :fascent, bbin = true)


    divergences[i] = div
    variances[i] = var

    println(i)
end


divergences = hcat(divergences...)
variances = hcat(variances...)

avg_divergence = sum(divergences, dims = 2) / run_number
avg_variance = sum(variances, dims = 2) / run_number



# ppm _ start _ spawns _ rate _ lineage_habit
name = "ppm_fbbin_10_25_fa_100"

using DelimitedFiles
writedlm("sub_bench_div_" * name * ".csv", divergences, ',')
writedlm("sub_bench_var_" * name * ".csv", variances, ",")



using Plots
plot([avg_divergence, avg_variance])