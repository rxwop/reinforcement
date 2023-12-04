cd = @__DIR__

using RxCiphers
text = Txt(read(raw"reproducible\\m10b.txt", String))
tokenise!(text)

include(joinpath(cd , "..", "spm reinforcement.jl"))
import Random.randperm, Statistics.mean


# Parameters of data collection
G_set = 600
run_number = 200

m = collect(1:130)
m = reshape(m, (5, :))
for (i, j) in enumerate([22, 18, 15, 14, 8])
    m[i, :] = circshift(m[i, :], - j)
end
m = m[[4,2,3,1,5], :]
target = invPermutation(vec(m))

inv_target = target.permutation

plain = target(text)
f_target = quadgramlog(plain)


print("Beginning runs...")
# 15.5%, 332.3 gen - :fascent w/ choice_weights
# 27.5%, 450.0 gen - :fascent w/ choice_weights, amnesia
# 97.5%, 372.93 gen - :ascent w/ choice_weights, amnesia
# 63.0%, 320.35 gen - :ascent w/ choice_weights
# 63.0%, 336.76 gen - :ascent
solve_log = Vector{Int}(undef, run_number)

for run in 1:run_number
    println("Run ", run)


    solver = PermutationSolve(
        text,
        130,
        150, # spawns
        10.; # learning rate (M)
        lineage_habit = :ascent,
        fscore = true
    )

    # initialise solve checker
    solved = false
    solve_num = G_set + 1

    for i in 1:G_set
        nextgen!(solver)

        # if lineage_certainty(solver.spm) > 0.95
        #     forget!(solver.spm, 0.8)
        # end
        
        if !solved
            if solver.parent.permutation == inv_target
                solve_num = i
                solved = true
            end
        end
    end


    solve_log[run] = solve_num

end

solve_percentage = count(!=(G_set + 1), solve_log) / run_number * 100

filter!(!=(G_set + 1), solve_log)

############# EXPORT DATA ###############



println(round(solve_percentage; digits = 1), "% of runs ended in solution")


if !isempty(solve_log)
    println("On average, solved in: ", mean(solve_log), " generations")
end

