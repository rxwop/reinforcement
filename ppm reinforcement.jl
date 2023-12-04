using RxCiphers
import LinearAlgebra.checksquare, Base.getindex, Base.setindex!
import StatsBase.sample, StatsBase.pweights


mutable struct PosProbMatrix{T <: AbstractFloat}
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

    function PosProbMatrix{T}(arr::Matrix{T}, size::Int, activation::Function, update_rate::T) where T <: AbstractFloat
        x = new{T}(arr)
        x.prob = normalise!(activation.(arr), 2)
        x.activation = activation
        x.needs_recalc = false

        x.size = size
        x.update_rate = update_rate
        x.indices = CartesianIndices((size, size))
        return x
    end
end

default_activation(f) = 1 + tanh(1000 * f)

function PosProbMatrix{T}(n::Int, update_rate::T; activation::Function = default_activation) ::PosProbMatrix{T} where T <: AbstractFloat
    return PosProbMatrix{T}(zeros(T, n, n), n, activation, update_rate)
end

function PosProbMatrix{T}(a::Matrix{T}, update_rate::T; activation::Function = default_activation) ::PosProbMatrix{T} where T <: AbstractFloat
    return PosProbMatrix{T}(a, checksquare(a), activation, update_rate)
end

function set_xy!(ppm::PosProbMatrix, x_to_y::Vector{Int})
    try ppm.xy
        println("WARNING: PosProbMatrix xy/yx mapping is already initialised, overwrite will occur.")
    catch err
    end

    ppm.xy = x_to_y
    ppm.yx = invperm(x_to_y)

    return ppm
end

function set_yx!(ppm::PosProbMatrix, y_to_x::Vector{Int})
    try ppm.yx
        println("WARNING: PosProbMatrix xy/yx mapping is already initialised, overwrite will occur.")
    catch err
    end

    ppm.yx = y_to_x
    ppm.xy = invperm(y_to_x)

    return ppm
end


function recalculate_prob!(ppm::PosProbMatrix) ::PosProbMatrix
    if ppm.needs_recalc
        ppm.prob = normalise!(ppm.activation.(ppm.arr), 2)
        ppm.needs_recalc = false
    end

    return ppm
end

getprob!(ppm::PosProbMatrix) = recalculate_prob!(ppm).prob
getprob(ppm::PosProbMatrix) = copy(getprob!(ppm))
getindex(ppm::PosProbMatrix, args...) = getindex(getprob!(ppm), args...)

function update!(ppm::PosProbMatrix{T}, x1::Int, x2::Int, y1::Int, y2::Int, dF::T) ::PosProbMatrix{T} where T <: AbstractFloat
    d1 = dF * ppm.update_rate
    d2 = dF * ppm.update_rate

    ppm.arr[x1, y1] -= d1
    ppm.arr[x1, y2] += d1

    ppm.arr[x2, y2] -= d2
    ppm.arr[x2, y1] += d2

    ppm.needs_recalc = true

    return ppm
end


# choice weights / bias
function draw(ppm::PosProbMatrix, number::Int, f::Function = identity) ::Vector{Tuple{Int, Int, Int, Int}}

    draw_matrix = f(getprob(ppm))

    for (x,y) in enumerate(ppm.xy)
        draw_matrix[x,y] = 0.
    end

    samples = sample(ppm.indices, pweights( draw_matrix ), number, replace = false) # Vector{Tuple} (x1, y2)

    x1 = [i[1] for i in samples]
    y2 = [i[2] for i in samples]

    x2 = [ppm.yx[i] for i in y2]
    y1 = [ppm.xy[i] for i in x1]


    return collect(zip(x1, x2, y1, y2))
end


function swap!(ppm::PosProbMatrix, (x1, x2, y1, y2)::Tuple{Int, Int, Int, Int})
    switch!(ppm.xy, x1, x2)
    switch!(ppm.yx, y1, y2)
    return ppm
end


function next_parent_dF(ppm::PosProbMatrix, swaps::Vector{Tuple{Int, Int, Int, Int}}, children::Vector{T}, child_dF::Vector{Float64}, lineage_habit = :ascent; parent::T) where T
    if !(length(swaps) == length(children) == length(child_dF))
        error("Each swap must have an associated child and delta fitness")
    end

    iter = zip(swaps, children, child_dF)

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
        swap!(ppm, best[1])
    end

    return best[2:3]
end



# Fitness function for PosProbMat, calculates Kullback-Leibler Divergence between PosProbMat and the target (certain) PosProbMat
function ppm_kldiv(ppm::PosProbMatrix, target_x_to_y::Vector{Int}) ::Float64
    s = 0

    for (x, y) in enumerate(target_x_to_y)
        s += log2(ppm[x, y])
    end

    return - s / ppm.size
end


# Fitness function for PosProbMat, calculates Kullback-Leibler Divergence between PosProbMat and the target (certain) PosProbMat
function ppm_kldiv(ppm::PosProbMatrix, target_x_to_y::Vector{Int}, ref_freq::Vector{Float64}) ::Float64
    s = 0

    for ((x, y), f) in zip(enumerate(target_x_to_y), ref_freq)
        s += f * log2(ppm[x, y])
    end

    return - s
end


function lineage_certainty(ppm::PosProbMatrix)
    avg = 0.

    for x in 1:ppm.size
        avg += ppm[x, ppm.xy[x]]
    end

    return avg / ppm.size
end




#################################################################################################################
















using Plots

# x (token) -> y (position in sub mapping)
# this solves for the inverse substitution, then inverts it
# Substitution solver, where ppM supervises a single substitution lineage

function substitution_solve(
    txt::Txt,
    generations::Int = 75,
    spawns::Int = 7,
    reinforce_rate::Float64 = 25.;
    lineage_habit::Symbol = :fascent,
    ref_freq::Vector{Float64} = monogram_freq
)
    solver = SubstitutionSolve(
        txt,
        frequency_matched_Substitution(txt, ref_freq),
        spawns,
        reinforce_rate;
        lineage_habit = lineage_habit,
        ref_freq = ref_freq,
        bbin = true
    )

    for i in 1:generations
        nextgen!(solver)
    end

    return invert!(solver.parent), solver.ppm
    # final Substitution in lineage // Reinforced Matrix
end












using LinearAlgebra
perm_matrix(permutation) = Matrix{Float64}(I, length(permutation), length(permutation))[invperm(permutation), :]
perm_matrix(perma, permb) = perm_matrix(invperm(perma)[checkperm(permb)])

mutable struct SubstitutionSolve
    text::Txt
    spawns::Int
    lineage_habit::Symbol

    ppm::PosProbMatrix

    parent::Substitution
    children::Vector{Substitution}

    fitness::Float64

    function SubstitutionSolve(
        text::Txt, start::Substitution, spawns::Int, reinforce_rate::Float64 = 0.5;
        lineage_habit::Symbol = :ascent,
        ref_freq::Vector{Float64} = monogram_freq,
        bbin::Bool = true
        )

        if bbin
            ppm = PosProbMatrix{Float64}(bbin_probabilities(text, ref_freq), reinforce_rate)
        else
            ppm = PosProbMatrix{Float64}(26, reinforce_rate)
        end
        set_yx!(ppm, start.mapping)

        new(
            text,
            spawns,
            lineage_habit,
            ppm,
            start,
            Vector{Substitution}(),
            quadgramlog(apply(start, text))
        )
    end
end

function nextgen!(s::SubstitutionSolve)

    swaps = draw(s.ppm, s.spawns)
    s.children = [switch(s.parent, y1, y2) for (x1, x2, y1, y2) in swaps]
    delta_F = [quadgramlog(new_sub(s.text)) for new_sub in s.children] .- s.fitness
    # generates new swaps from ppM and calculates dF


    for ((x1, x2, y1, y2), dF) in zip(swaps, delta_F) # Update P with ALL the data
        update!(s.ppm, x1, x2, y1, y2, dF)
    end


    s.parent, dF = next_parent_dF(s.ppm, swaps, s.children, delta_F, s.lineage_habit; parent = s.parent)
    s.fitness += dF
    # advance lineage
end