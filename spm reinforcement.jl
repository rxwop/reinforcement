using RxCiphers
import LinearAlgebra.checksquare, Base.getindex, Base.setindex!
import StatsBase.sample, StatsBase.pweights
import Statistics.mean


function normalise_nodiag!(arr::Matrix{AbstractFloat}, dims) ::Matrix{AbstractFloat}
    arr ./= (sum(arr; dims = dims) - sum(diag(arr)))
    return arr
end

mutable struct SsrProbMatrix{T <: AbstractFloat}
    arr::Matrix{T}

    # activation must map : [-Inf, Inf] -> [0, Inf]
    activation::Function
    prob::Matrix{T}
    needs_recalc::Bool
    
    xy::Vector{Int}
    yx::Vector{Int}

    size::Int
    update_rate::T

    indices::CartesianIndices{2, Tuple{Base.OneTo{Int64}, Base.OneTo{Int64}}}

    function SsrProbMatrix{T}(arr::Matrix{Float64}, size::Int, activation::Function, update_rate::Float64)
        x = new{T}(arr)
        x.prob = normalise_nodiag!(activation.(arr), 2)
        for i in 1:size
            x.prob[i,i] = 0.
        end

        x.activation = activation
        x.needs_recalc = false

        x.size = size
        x.update_rate = update_rate
        x.indices = CartesianIndices((size, size))
        return x
    end
end

default_activation(f) = 1 + tanh(1000 * f)

function SsrProbMatrix{T}(n::Int, update_rate::T; activation::Function = default_activation) ::SsrProbMatrix{T} where T <: AbstractFloat
    return SsrProbMatrix{T}(zeros(T, n, n), n, activation, update_rate)
end

function SsrProbMatrix{T}(a::Matrix{Float64}, update_rate::Float64; activation::Function = default_activation) ::SsrProbMatrix{T} where T <: AbstractFloat
    return SsrProbMatrix{T}(a, checksquare(a), activation, update_rate)
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

function recalculate_prob!(spm::SsrProbMatrix) ::SsrProbMatrix
    if spm.needs_recalc
        spm.prob = normalise_nodiag!(spm.activation.(spm.arr), 2)
        spm.needs_recalc = false
    end

    return spm
end

getprob!(spm::SsrProbMatrix) = recalculate_prob!(spm).prob
getprob(spm::SsrProbMatrix) = copy(getprob!(spm))
getindex(spm::SsrProbMatrix, args...) = getindex(getprob!(spm), args...)


function update!(spm::SsrProbMatrix, x1::Int, x2::Int, x3::Int, y1::Int, y2::Int, y3::Int, dF::Float64) ::SsrProbMatrix
    d1 = dF * spm.update_rate
    d2 = dF * spm.update_rate
    d3 = dF * spm.update_rate

    spm.arr[x1, y1] -= d1
    spm.arr[x1, y2] += d1

    spm.arr[x2, y2] -= d2
    spm.arr[x2, y1] += d2

    spm.arr[x3, y3] -= d3
    spm.arr[x3, y1] += d3

    spm.needs_recalc = true

    return spm
end


bernoulli(w::Float64) = rand() < w
import Statistics.stdm
function draw(spm::SsrProbMatrix, number::Int, f::Function = identity) ::Vector{NTuple{6, Int}}
    
    draw_matrix = f(getprob(spm))

    for (x,y) in enumerate(spm.xy)
        draw_matrix[x,y] = 0.
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

function lineage_certainty(spm::SsrProbMatrix)
    avg = 0.

    for x in 1:spm.size
        avg += spm[x, spm.xy[x]]
    end

    return avg / spm.size
end


#################################################################################################################


function forget!(spm::SsrProbMatrix, a::Float64)
    spm.arr *= a
    return spm
end









using Plots

# x (token) -> y (successor in permutation)
# this solves for the inverse permutation, then inverts it

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



mutable struct PermutationSolve
    text::Txt
    spawns::Int
    lineage_habit::Symbol
    N::Int

    amnesia_coeff::Float64
    amnesia_threshold::Float64
    has_amnesia::Bool

    spm::SsrProbMatrix

    parent::Permutation
    children::Vector{Permutation}

    fitness::Float64

    function PermutationSolve(
        text::Txt, N::Int, spawns::Int, reinforce_rate::Float64 = 0.5;
        lineage_habit::Symbol = :ascent,
        fscore::Bool = true,
        amnesia_threshold::Float64 = 0.95,
        amnesia_coeff::Float64 = 1.0
        )

        e = ArgumentError("Amnesia coefficient must be positive and not greater than 1.0")
        if !(0.0 < amnesia_coeff <= 1.0)
            throw(e)
        elseif amnesia_coeff == 1.0
            has_amnesia = false
        elseif amnesia_threshold >= 1.0
            has_amnesia = false
        else
            has_amnesia = true
        end

        S = N + 1

        start = invPermutation(N, true) # guesses INVERSE permutation

        if fscore
            follow_scr = exp10.(follow_scores(RxCiphers.unshape(text.tokenised, start)))
            blanks = ones(Float64, (S, S)) # For the extra split token
            blanks[1:N, 1:N] = follow_scr # stitch in the follow matrix
            blanks[1:N, S] = mean(follow_scr; dims = 2)
            spm = SsrProbMatrix{Float64}(blanks, reinforce_rate)
            
        else
            spm = SsrProbMatrix{Float64}(S, reinforce_rate)
        end
        set_xy!(spm, splitsuccessors(start.permutation))

        new(
            text,
            spawns,
            lineage_habit,
            S,
            amnesia_coeff,
            amnesia_threshold,
            has_amnesia,
            spm,
            start,
            Vector{Permutation}(),
            quadgramlog(apply(start, text))
        )
    end
end


function nextgen!(s::PermutationSolve)
    splices = draw(s.spm, s.spawns)

    parent_split_perm = [s.N ; s.parent.permutation]
    s.children = [invPermutation(RxCiphers.circorigin(RxCiphers.splice(parent_split_perm, 
    findfirst(==(y1), parent_split_perm), 
    findfirst(==(y2), parent_split_perm), 
    findfirst(==(x3), parent_split_perm)), s.N)[2:end], true) for (x1, x2, x3, y1, y2, y3) in splices]

    delta_F = [quadgramlog(new_p(s.text)) for new_p in s.children] .- s.fitness
    # generates new splices from ppM and calculates dF


    for ((x1, x2, x3, y1, y2, y3), dF) in zip(splices, delta_F) # Update P with ALL the data
        update!(s.spm, x1, x2, x3, y1, y2, y3, dF)
    end


    s.parent, dF = next_parent_dF(s.spm, splices, s.children, delta_F, s.lineage_habit; parent = s.parent)
    s.fitness += dF
    # advance lineage

    if s.has_amnesia
        if lineage_certainty(s.spm) > s.amnesia_threshold
            forget!(s.spm, s.amnesia_coeff)
        end
    end
end




















mutable struct MultPermutationSolve
    lineage_num::Int

    text::Txt
    spawns::Int
    lineage_habit::Symbol
    N::Int

    spm::SsrProbMatrix

    parent::Vector{Permutation}
    children::Matrix{Permutation}

    fitness::Float64

    function MultPermutationSolve(
        lineage_num::Int, text::Txt, N::Int, spawns::Int, reinforce_rate::Float64 = 0.5;
        lineage_habit::Symbol = :ascent,
        fscore::Bool = true
        )

        S = N + 1

        starts = [invPermutation(N, true)] # guesses INVERSE permutation

        if fscore
            follow_scr = exp10.(follow_scores(RxCiphers.unshape(text.tokenised, start)))
            blanks = fill(1e-10, (S, S)) # For the extra split token
            blanks[1:N, 1:N] = follow_scr # stitch in the follow matrix
            spm = SsrProbMatrix(blanks, reinforce_rate)
            
        else
            spm = SsrProbMatrix(S, reinforce_rate)
        end
        set_xy!(spm, splitsuccessors(start.permutation))

        new(
            text,
            spawns,
            lineage_habit,
            S,
            spm,
            starts,
            Matrix{Permutation}(),
            quadgramlog(apply(start, text))
        )
    end
end

function nextgen!(s::MultPermutationSolve)

    splices = draw(s.spm, s.spawns)

    parent_split_perm = [s.N ; s.parent.permutation]
    new_perms = [invPermutation(RxCiphers.circorigin(RxCiphers.splice(parent_split_perm, 
    findfirst(==(y1), parent_split_perm), 
    findfirst(==(y2), parent_split_perm), 
    findfirst(==(x3), parent_split_perm)), s.N)[2:end], true) for (x1, x2, x3, y1, y2, y3) in splices]

    delta_F = [quadgramlog(new_p(s.text)) for new_p in new_perms] .- s.fitness
    # generates new splices from ppM and calculates dF


    for ((x1, x2, x3, y1, y2, y3), dF) in zip(splices, delta_F) # Update P with ALL the data
        update!(s.spm, x1, x2, x3, y1, y2, y3, dF)
    end


    s.parent, dF = next_parent_dF(s.spm, splices, new_perms, delta_F, s.lineage_habit; parent = s.parent)
    s.fitness += dF
    # advance lineage
end









using LinearAlgebra
perm_matrix(permutation) = Matrix{Float64}(I, length(permutation), length(permutation))[checkperm(permutation), :]
perm_matrix(perma, permb) = perm_matrix(invperm(perma)[checkperm(permb)])