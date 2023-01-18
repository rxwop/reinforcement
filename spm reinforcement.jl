include("permutation.jl")
include("tuco.jl")
import LinearAlgebra.checksquare, Base.getindex, Base.setindex!
import StatsBase.sample, StatsBase.pweights


mutable struct SsrProbMatrix
    arr::Matrix{Float64}
    xy::Vector{Int}
    yx::Vector{Int}

    size::Int
    update_rate::Float64

    indices::CartesianIndices{2, Tuple{Base.OneTo{Int64}, Base.OneTo{Int64}}}

    function SsrProbMatrix(arr::Matrix{Float64}, size::Int, update_rate::Float64)
        x = new(arr)
        x.size = size
        x.update_rate = update_rate
        x.indices = CartesianIndices((size, size))
        return x
    end
end

function SsrProbMatrix(n::Int, update_rate::Float64) ::SsrProbMatrix
    return SsrProbMatrix(fill(1 / n, (n,n)), n, update_rate)
end

function SsrProbMatrix(a::Matrix{Float64}, update_rate::Float64) ::SsrProbMatrix
    return SsrProbMatrix(normalise!(clamp.(a, 0., 1.); dims = 2), checksquare(a), update_rate)
end

function set_xy!(spm::SsrProbMatrix, x_to_y::Vector{Int})
    try spm.xy
        println("WARNING: SsrProbMatrix xy/yx mapping is already initialised, overwrite will occur.")
    catch err
    end

    spm.xy = x_to_y
    spm.yx = invperm(x_to_y)

    return spm
end

function set_yx!(spm::SsrProbMatrix, y_to_x::Vector{Int})
    try spm.yx
        println("WARNING: SsrProbMatrix xy/yx mapping is already initialised, overwrite will occur.")
    catch err
    end

    spm.yx = y_to_x
    spm.xy = invperm(y_to_x)

    return spm
end




getindex(spm::SsrProbMatrix, args...) = getindex(spm.arr, args...)


function update!(spm::SsrProbMatrix, x1::Int, x2::Int, x3::Int, y1::Int, y2::Int, y3::Int, dF::Float64) ::SsrProbMatrix
    d1 = update_delta(dF, spm[x1, y1], spm[x1, y2], spm.update_rate)
    d2 = update_delta(dF, spm[x2, y2], spm[x2, y3], spm.update_rate)
    d3 = update_delta(dF, spm[x3, y3], spm[x3, y1], spm.update_rate)

    spm.arr[x1, y1] -= d1
    spm.arr[x1, y2] += d1

    spm.arr[x2, y2] -= d2
    spm.arr[x2, y1] += d2

    spm.arr[x3, y3] -= d3
    spm.arr[x3, y1] += d3

    for i in (x1, x2, x3)
        for j in (y1, y2, y3)
            spm.arr[i,j] = clamp(spm.arr[i,j], 0., 1.)
        end
    end

    return spm
end

# ZEROS ARE STUCK
# Limiting function taking (-inf, inf) -> (-p_new, p_old), with update_delta(0) == 0
function update_delta(delta_fitness::Float64, p_old::Float64, p_new::Float64, rate::Float64) ::Float64

    if p_old == 0. || p_new == 0.
        return 0.
    end

    if rate == Inf # maximum learn rate
        if delta_fitness == 0.
            return 0.0
        elseif delta_fitness > 0.
            return p_old
        else
            return (- p_new)
        end
    end


    mu = p_old + p_new # division by 2 carried out at the end
    diff = p_old - p_new # division by 2 carried out at the end

    z = coth(rate * delta_fitness)

    delta = 2 * p_old * p_new / (mu * z - diff)

    return clamp(delta, - p_new, p_old)
end


bernoulli(w::Float64) = rand() < w

function draw(spm::SsrProbMatrix, number::Int, f::Function = identity) ::Vector{NTuple{6, Int}}
    
    draw_matrix = f(copy(spm.arr))

    for (x,y) in enumerate(spm.xy)
        draw_matrix[x,y] = 0.
    end

    for i in 1:spm.size
        draw_matrix[i, i] = 0. # dodgy
    end

    samples = sample(spm.indices, pweights( draw_matrix ), number, replace = false) # Vector{Tuple} (x1, y2)

    x1 = [i[1] for i in samples]
    y2 = [i[2] for i in samples]

    x2 = [spm.yx[i] for i in y2]
    y1 = [spm.xy[i] for i in x1]

    x3 = copy(y2)
    y3 = [spm.xy[i] for i in x3]

    for i in 1:number # for each x3, y3
        lim = x1[i]

        while bernoulli(spm[x3[i], y3[i]]) # while bernoulli success for P[item and successor]
            if y3[i] == lim
                break
            else
                x3[i] = y3[i] # advance both x3, y3 in succession order
                y3[i] = spm.xy[y3[i]]
            end
        end
    end

    return collect(zip(x1, x2, x3, y1, y2, y3))
end


function next_parent_dF(spm::SsrProbMatrix, splices::Vector{NTuple{6, Int}}, children::Vector{T}, child_dF::Vector{Float64}, lineage_habit = :ascent; parent::T) where T
    if !(length(splices) == length(children) == length(child_dF))
        error("Each swap must have an associated child and delta fitness")
    end

    iter = zip(splices, children, child_dF)

    # Set new parent and fitness
    if lineage_habit == :ascent
        best = first(iter)

        for item in iter # choose best child
            if item[3] > best[3]
                best = item
            end
        end

    elseif lineage_habit == :fascent
        best = (nothing, parent, 0.)

        for item in iter # choose best child
            if item[3] > best[3]
                best = item
            end
        end

    elseif lineage_habit == :random
        best = first(iter)

    elseif lineage_habit == :descent
        best = first(iter)

        for item in iter # choose worst child
            if item[3] < best[3]
                best = item
            end
        end

    elseif lineage_habit == :stationary
        best = (nothing, parent, 0.)

    else
        error("Invalid lineage habit kwarg")
    end

    if !isnothing(best[1])
        splice!(spm, best[1])
    end

    return best[2:3]
end



# Fitness function for SsrProbMat, calculates Kullback-Leibler Divergence between SsrProbMat and the target (certain) SsrProbMat
function spm_kldiv(spm::SsrProbMatrix, target_x_to_y::Vector{Int}) ::Float64
    s = 0

    for (x, y) in enumerate(target_x_to_y)
        s += log2(spm[x, y])
    end

    return - s / spm.size
end


# Fitness function for SsrProbMat, calculates Kullback-Leibler Divergence between SsrProbMat and the target (certain) SsrProbMat
function spm_kldiv(spm::SsrProbMatrix, target_x_to_y::Vector{Int}, ref_freq::Vector{Float64}) ::Float64
    s = 0

    for ((x, y), f) in zip(enumerate(target_x_to_y), ref_freq)
        s += f * log2(spm[x, y])
    end

    return - s
end


###########################################################################################################


function switch_231!(self::AbstractVector, posa::Int, posb::Int, posc::Int) ::AbstractVector
    self[posa], self[posb], self[posc] = self[posb], self[posc], self[posa]
    return self
end
switch_231(v::Vector{Int}, posa::Int, posb::Int, posc::Int) ::Vector{Int} = switch_231!(copy(v), posa, posb, posc)

function switch_312!(self::AbstractVector, posa::Int, posb::Int, posc::Int) ::AbstractVector
    self[posa], self[posb], self[posc] = self[posc], self[posa], self[posb]
    return self
end
switch_312(v::Vector{Int}, posa::Int, posb::Int, posc::Int) ::Vector{Int} = switch_312!(copy(v), posa, posb, posc)


function splice!(spm::SsrProbMatrix, (x1, x2, x3, y1, y2, y3)::Tuple{Int, Int, Int, Int, Int, Int})
    switch_231!(spm.xy, x1, x2, x3)
    switch_312!(spm.yx, y1, y2, y3)
    return spm
end



function successors(permutation::Vector{Int}) ::Vector{Int}
    L = length(checkperm(permutation))
    out = Vector{Int}(undef, L)

    for i in 1:L-1
        out[permutation[i]] = permutation[i + 1]
    end

    out[permutation[L]] = permutation[1]

    return out
end

function splitsuccessors(permutation::Vector{Int}) ::Vector{Int}
    L = length(checkperm(permutation))
    split = L + 1
    out = Vector{Int}(undef, split)

    for i in 1:L-1
        out[permutation[i]] = permutation[i + 1]
    end


    out[permutation[L]] = split
    out[split] = permutation[1]

    return out
end

###########################################################################################################

function bigram_follow_score(cola::Vector{Int}, colb::Vector{Int}) ::Float64
    s = 0.
    n = 0

    for bigram in eachrow([cola colb])
        if any(bigram .== 0)
            continue
        end

        s += bigram_scores[bigram[1], bigram[2]]
        n += 1
    end

    return s / n
end

function follow_scores(matrix::Matrix{Int}) ::Matrix{Float64}
    n = size(matrix)[1]
    follow_scores = Matrix{Float64}(undef, (n, n))

    for i in 1:n
        for j in 1:n
            follow_scores[i, j] = bigram_follow_score(matrix[i, :], matrix[j, :])
        end
    end

    return follow_scores
end

#################################################################################################################
















using Plots

# x (token) -> y (successor in permutation)
# this solves for the inverse permutation, then inverts it
function debug_permutation_solve(
    target::Permutation,
    txt::Txt,
    generations::Int,
    spawns::Int,
    reinforce_rate::Float64 = 0.5;
    lineage_habit::Symbol = :ascent,
    frame_skip::Int = 1,
)
    N = target.n
    S = N + 1
    target = invert(target)
    fitness_log = Vector{Float64}(undef, generations + 1)
    div_log = Vector{Float64}(undef, generations + 1)
    target_xy = splitsuccessors(target.permutation)



    parent_perm = invPermutation(N, true) # guesses INVERSE permutation

    follow_scr = exp10.(follow_scores(unshape(txt.tokenised, parent_perm)))
    blanks = fill(1e-10, (S, S))
    blanks[1:N, 1:N] = follow_scr
    spm = SsrProbMatrix(blanks, reinforce_rate)
    set_xy!(spm, splitsuccessors(parent_perm.permutation))



    parent_fitness = quadgramlog(apply(parent_perm, txt))
    fitness_log[1] = parent_fitness
    div_log[1] = spm_kldiv(spm, target_xy)

    gr(format=:png)
    p = heatmap(spm.arr, clims = (0, 1), aspect_ratio = :equal, xlabel = "x", ylabel = "y", xticks = false, yticks = false)

    anim = @animate for gen in 1:generations
        println(gen)
        splices = draw(spm, spawns)

        parent_split_perm = [S ; parent_perm.permutation]
        new_perms = [invPermutation(circorigin(splice(parent_split_perm, 
        findfirst(==(y1), parent_split_perm), 
        findfirst(==(y2), parent_split_perm), 
        findfirst(==(x3), parent_split_perm)), S)[2:end], true) for (x1, x2, x3, y1, y2, y3) in splices]

        delta_F = [quadgramlog(new_p(txt)) for new_p in new_perms] .- parent_fitness
        # generates new splices from ppM and calculates dF


        for ((x1, x2, x3, y1, y2, y3), dF) in zip(splices, delta_F) # Update P with ALL the data
            update!(spm, x1, x2, x3, y1, y2, y3, dF)
        end


        parent_perm, dF = next_parent_dF(spm, splices, new_perms, delta_F, lineage_habit; parent = parent_perm)
        parent_fitness += dF
        # advance lineage


        fitness_log[gen + 1] = parent_fitness
        div_log[gen + 1] = spm_kldiv(spm, target_xy)

        p[1][1][:z] = spm.arr

    end every frame_skip

    
    gif(anim, "anim.gif")

    invert!(parent_perm) # return the "solved" FORWARDS permutation

    return spm, parent_perm, fitness_log, div_log
    # Reinforced Matrix // final Substitution in lineage // Vector of fitness values vs generations // Vector of divergence vs gen
end


function permutation_solve(
    N::Int,
    txt::Txt,
    generations::Int,
    spawns::Int,
    reinforce_rate::Float64 = 0.5;
    lineage_habit::Symbol = :ascent,
)
    S = N + 1

    parent_perm = invPermutation(N, true) # guesses INVERSE permutation

    follow_scr = exp10.(follow_scores(unshape(txt.tokenised, parent_perm)))
    blanks = fill(1e-10, (S, S))
    blanks[1:N, 1:N] = follow_scr
    spm = SsrProbMatrix(blanks, reinforce_rate)
    set_xy!(spm, splitsuccessors(parent_perm.permutation))



    parent_fitness = quadgramlog(apply(parent_perm, txt))

    for gen in 1:generations
        splices = draw(spm, spawns)

        parent_split_perm = [S ; parent_perm.permutation]
        new_perms = [invPermutation(circorigin(splice(parent_split_perm, 
        findfirst(==(y1), parent_split_perm), 
        findfirst(==(y2), parent_split_perm), 
        findfirst(==(x3), parent_split_perm)), S)[2:end], true) for (x1, x2, x3, y1, y2, y3) in splices]

        delta_F = [quadgramlog(new_p(txt)) for new_p in new_perms] .- parent_fitness
        # generates new splices from ppM and calculates dF


        for ((x1, x2, x3, y1, y2, y3), dF) in zip(splices, delta_F) # Update P with ALL the data
            update!(spm, x1, x2, x3, y1, y2, y3, dF)
        end


        parent_perm, dF = next_parent_dF(spm, splices, new_perms, delta_F, lineage_habit; parent = parent_perm)
        parent_fitness += dF
        # advance lineage

    end

    invert!(parent_perm) # return the "solved" FORWARDS permutation

    return spm, parent_perm
    # Reinforced Matrix // final Substitution in lineage
end











using LinearAlgebra
perm_matrix(permutation) = Matrix{Float64}(I, length(permutation), length(permutation))[checkperm(permutation), :]

function permutation_show_off(
    target::Permutation,
    txt::Txt,
    generations::Int,
    spawns::Int,
    reinforce_rate::Float64 = 0.5;
    lineage_habit::Symbol = :ascent,
    frame_skip::Int = 1,
)
    N = target.n
    S = N + 1
    target = invert(target)
    fitness_log = Vector{Float64}(undef, generations + 1)
    div_log = Vector{Float64}(undef, generations + 1)
    target_xy = splitsuccessors(target.permutation)



    parent_perm = invPermutation(N, true) # guesses INVERSE permutation

    follow_scr = exp10.(follow_scores(unshape(txt.tokenised, parent_perm)))
    blanks = fill(1e-10, (S, S))
    blanks[1:N, 1:N] = follow_scr
    spm = SsrProbMatrix(blanks, reinforce_rate)
    set_xy!(spm, splitsuccessors(parent_perm.permutation))



    parent_fitness = quadgramlog(apply(parent_perm, txt))
    fitness_log[1] = parent_fitness
    div_log[1] = spm_kldiv(spm, target_xy)

    gr(format=:png)
    p1 = heatmap(spm.arr, clims = (0, 1), c = cgrad(["0x131314", "0x0048ff"]), aspect_ratio = :equal, xlabel = "x", ylabel = "y", xticks = false, yticks = false, title = "spM")
    p2 = heatmap(perm_matrix(parent_perm.permutation), clims = (0, 1), c = cgrad(["0x131314", "0x30b891"]), legend = :none, aspect_ratio = :equal, xlabel = "x", ylabel = "y", xticks = false, yticks = false, title = "Permutation Lineage")
    p = plot(p1, p2)

    anim = @animate for gen in 1:generations
        println(gen)
        splices = draw(spm, spawns)

        parent_split_perm = [S ; parent_perm.permutation]
        new_perms = [invPermutation(circorigin(splice(parent_split_perm, 
        findfirst(==(y1), parent_split_perm), 
        findfirst(==(y2), parent_split_perm), 
        findfirst(==(x3), parent_split_perm)), S)[2:end], true) for (x1, x2, x3, y1, y2, y3) in splices]

        delta_F = [quadgramlog(new_p(txt)) for new_p in new_perms] .- parent_fitness
        # generates new splices from ppM and calculates dF


        for ((x1, x2, x3, y1, y2, y3), dF) in zip(splices, delta_F) # Update P with ALL the data
            update!(spm, x1, x2, x3, y1, y2, y3, dF)
        end


        parent_perm, dF = next_parent_dF(spm, splices, new_perms, delta_F, lineage_habit; parent = parent_perm)
        parent_fitness += dF
        # advance lineage


        fitness_log[gen + 1] = parent_fitness
        div_log[gen + 1] = spm_kldiv(spm, target_xy)

        p[1][1][:z] = spm.arr
        p[2][1][:z] = perm_matrix(parent_perm.permutation)

    end every frame_skip

    
    gif(anim, "anim.gif")

    invert!(parent_perm) # return the "solved" FORWARDS permutation

    return spm, parent_perm, fitness_log, div_log
    # Reinforced Matrix // final Substitution in lineage // Vector of fitness values vs generations // Vector of divergence vs gen
end