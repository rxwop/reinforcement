include("substitution.jl")
include("tuco.jl")
import LinearAlgebra.checksquare, Base.getindex, Base.setindex!
import StatsBase.sample, StatsBase.pweights


mutable struct PosProbMatrix
    arr::Matrix{Float64}
    xy::Vector{Int}
    yx::Vector{Int}

    size::Int
    update_rate::Float64

    indices::CartesianIndices{2, Tuple{Base.OneTo{Int64}, Base.OneTo{Int64}}}

    function PosProbMatrix(arr::Matrix{Float64}, size::Int, update_rate::Float64)
        x = new(arr)
        x.size = size
        x.update_rate = update_rate
        x.indices = CartesianIndices((size, size))
        return x
    end
end

function PosProbMatrix(n::Int, update_rate::Float64) ::PosProbMatrix
    return PosProbMatrix(fill(1 / n, (n,n)), n, update_rate)
end

function PosProbMatrix(a::Matrix{Float64}, update_rate::Float64) ::PosProbMatrix
    return PosProbMatrix(normalise!(clamp.(a, 0., 1.); dims = 2), checksquare(a), update_rate)
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




getindex(ppm::PosProbMatrix, args...) = getindex(ppm.arr, args...)


function update!(ppm::PosProbMatrix, x1::Int, x2::Int, y1::Int, y2::Int, dF::Float64) ::PosProbMatrix
    d1 = update_delta(dF, ppm[x1, y1], ppm[x1, y2], ppm.update_rate)
    d2 = update_delta(dF, ppm[x2, y2], ppm[x2, y1], ppm.update_rate)

    ppm.arr[x1, y1] -= d1
    ppm.arr[x1, y2] += d1

    ppm.arr[x2, y2] -= d2
    ppm.arr[x2, y1] += d2

    for i in (x1, x2)
        for j in (y1, y2)
            ppm.arr[i,j] = clamp(ppm.arr[i,j], 0., 1.)
        end
    end

    return ppm
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


# choice weights / bias
function draw(ppm::PosProbMatrix, number::Int, f::Function = identity) ::Vector{Tuple{Int, Int, Int, Int}}

    draw_matrix = f(copy(ppm.arr))

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






#################################################################################################################
















using Plots

# x (token) -> y (position in sub mapping)
# this solves for the inverse substitution, then inverts it
# Substitution solver, where ppM supervises a single substitution lineage
function debug_substitution_solve(
    target::Substitution,
    txt::Txt,
    generations::Int,
    spawns::Int,
    reinforce_rate::Float64 = 0.5;
    lineage_habit::Symbol = :ascent,
    frame_skip::Int = 1,
    ref_freq::Vector{Float64} = monogram_freq
)

    fitness_log = Vector{Float64}(undef, generations + 1)
    div_log = Vector{Float64}(undef, generations + 1)



    parent_sub = frequency_matched_Substitution(txt, ref_freq) # guesses FORWARDS substitution
    invert!(parent_sub)

    ppm = PosProbMatrix(bbin_probabilities(txt, ref_freq), reinforce_rate)
    set_yx!(ppm, parent_sub.mapping)



    parent_fitness = quadgramlog(apply(parent_sub, txt))
    fitness_log[1] = parent_fitness
    div_log[1] = ppm_kldiv(ppm, target.mapping)

    gr(format=:png)
    p = heatmap(ppm.arr, clims = (0, 1), aspect_ratio = :equal, xlabel = "x", ylabel = "y", xticks = false, yticks = false)

    anim = @animate for gen in 1:generations
        println(gen)
        swaps = draw(ppm, spawns)
        new_substitutions = [switch(parent_sub, y1, y2) for (x1, x2, y1, y2) in swaps]
        delta_F = [quadgramlog(new_sub(txt)) for new_sub in new_substitutions] .- parent_fitness
        # generates new swaps from ppM and calculates dF


        for ((x1, x2, y1, y2), dF) in zip(swaps, delta_F) # Update P with ALL the data
            update!(ppm, x1, x2, y1, y2, dF)
        end


        parent_sub, dF = next_parent_dF(ppm, swaps, new_substitutions, delta_F, lineage_habit; parent = parent_sub)
        parent_fitness += dF
        # advance lineage


        fitness_log[gen + 1] = parent_fitness
        div_log[gen + 1] = ppm_kldiv(ppm, target.mapping)

        p[1][1][:z] = ppm.arr
        # scatter!([i[4] for i in swaps], [i[1] for i in swaps], color = :blue, markershape = :rect, label = false)
        # scatter!([i[3] for i in swaps], [i[2] for i in swaps], color = :blue, markershape = :rect, label = false)
        # scatter!([i[3] for i in swaps], [i[1] for i in swaps], color = :blue, markershape = :rect, label = false)
        # scatter!([i[4] for i in swaps], [i[2] for i in swaps], color = :blue, markershape = :rect, label = false)
    end every frame_skip

    
    gif(anim, "anim.gif")

    invert!(parent_sub) # return the "solved" FORWARDS substitution

    return ppm, parent_sub, fitness_log, div_log
    # Reinforced Matrix // final Substitution in lineage // Vector of fitness values vs generations // Vector of divergence vs gen
end









function substitution_solve(
    txt::Txt,
    generations::Int,
    spawns::Int,
    reinforce_rate::Float64 = 0.5;
    lineage_habit::Symbol = :ascent,
    ref_freq::Vector{Float64} = monogram_freq
)


    parent_sub = frequency_matched_Substitution(txt, ref_freq) # guesses FORWARDS substitution
    invert!(parent_sub)

    ppm = PosProbMatrix(bbin_probabilities(txt, ref_freq), reinforce_rate)
    set_yx!(ppm, parent_sub.mapping)



    parent_fitness = quadgramlog(apply(parent_sub, txt))

    gr(format=:png)
    p = heatmap(ppm.arr, clims = (0, 1), aspect_ratio = :equal, xlabel = "x", ylabel = "y", xticks = false, yticks = false)

    for gen in 1:generations
        swaps = draw(ppm, spawns)
        new_substitutions = [switch(parent_sub, y1, y2) for (x1, x2, y1, y2) in swaps]
        delta_F = [quadgramlog(new_sub(txt)) for new_sub in new_substitutions] .- parent_fitness
        # generates new swaps from ppM and calculates dF


        for ((x1, x2, y1, y2), dF) in zip(swaps, delta_F) # Update P with ALL the data
            update!(ppm, x1, x2, y1, y2, dF)
        end


        parent_sub, dF = next_parent_dF(ppm, swaps, new_substitutions, delta_F, lineage_habit; parent = parent_sub)
        parent_fitness += dF
        # advance lineage
    end


    invert!(parent_sub) # return the "solved" FORWARDS substitution

    return ppm, parent_sub
    # Reinforced Matrix // final Substitution in lineage
end












using LinearAlgebra
perm_matrix(permutation) = Matrix{Float64}(I, length(permutation), length(permutation))[invperm(permutation), :]

function substitution_show_off(
    target::Substitution,
    txt::Txt,
    generations::Int,
    spawns::Int,
    reinforce_rate::Float64 = 0.5;
    lineage_habit::Symbol = :ascent,
    frame_skip::Int = 1,
    ref_freq::Vector{Float64} = monogram_freq
)

    fitness_log = Vector{Float64}(undef, generations + 1)
    div_log = Vector{Float64}(undef, generations + 1)



    parent_sub = Substitution(target.size) # guesses FORWARDS substitution

    ppm = PosProbMatrix(bbin_probabilities(txt, ref_freq), reinforce_rate)
    set_yx!(ppm, parent_sub.mapping)



    parent_fitness = quadgramlog(apply(parent_sub, txt))
    fitness_log[1] = parent_fitness
    div_log[1] = ppm_kldiv(ppm, target.mapping)

    gr(format=:png)
    p1 = heatmap(ppm.arr, clims = (0, 1), c = cgrad(["0x131314", "0xf02259"]), aspect_ratio = :equal, xlabel = "x", ylabel = "y", xticks = false, yticks = false, title = "ppM")
    p2 = heatmap(perm_matrix(parent_sub.mapping), clims = (0, 1.3), c = cgrad(["0x131314", "0xe6e667"]), legend = :none, aspect_ratio = :equal, xlabel = "x", ylabel = "y", xticks = false, yticks = false, title = "Substitution Lineage")
    p = plot(p1, p2)

    anim = @animate for gen in 1:generations
        println(gen)
        swaps = draw(ppm, spawns)
        new_substitutions = [switch(parent_sub, y1, y2) for (x1, x2, y1, y2) in swaps]
        delta_F = [quadgramlog(new_sub(txt)) for new_sub in new_substitutions] .- parent_fitness
        # generates new swaps from ppM and calculates dF


        for ((x1, x2, y1, y2), dF) in zip(swaps, delta_F) # Update P with ALL the data
            update!(ppm, x1, x2, y1, y2, dF)
        end


        parent_sub, dF = next_parent_dF(ppm, swaps, new_substitutions, delta_F, lineage_habit; parent = parent_sub)
        parent_fitness += dF
        # advance lineage


        fitness_log[gen + 1] = parent_fitness
        div_log[gen + 1] = ppm_kldiv(ppm, target.mapping)

        p[1][1][:z] = ppm.arr
        p[2][1][:z] = perm_matrix(parent_sub.mapping)

    end every frame_skip

    
    gif(anim, "anim.gif")

    invert!(parent_sub) # return the "solved" FORWARDS substitution

    return ppm, parent_sub, fitness_log, div_log
    # Reinforced Matrix // final Substitution in lineage // Vector of fitness values vs generations // Vector of divergence vs gen
end





function bench_substitution_solve(
    start::Substitution,
    inv_target::Substitution,
    txt::Txt,
    generations::Int,
    spawns::Int,
    reinforce_rate::Float64 = 0.5;
    lineage_habit::Symbol = :ascent,
    ref_freq::Vector{Float64} = monogram_freq,
    bbin::Bool = false
)

    fitness_log = Vector{Float64}(undef, generations + 1)
    solved = false
    solve_num = generations


    parent_sub = start # guesses FORWARDS substitution

    if bbin
        ppm = PosProbMatrix(bbin_probabilities(txt, ref_freq), reinforce_rate)
    else
        ppm = PosProbMatrix(26, reinforce_rate)
    end
    set_yx!(ppm, parent_sub.mapping)



    parent_fitness = quadgramlog(apply(parent_sub, txt))
    fitness_log[1] = parent_fitness


    for gen in 1:generations
        swaps = draw(ppm, spawns)
        new_substitutions = [switch(parent_sub, y1, y2) for (x1, x2, y1, y2) in swaps]
        delta_F = [quadgramlog(new_sub(txt)) for new_sub in new_substitutions] .- parent_fitness
        # generates new swaps from ppM and calculates dF


        for ((x1, x2, y1, y2), dF) in zip(swaps, delta_F) # Update P with ALL the data
            update!(ppm, x1, x2, y1, y2, dF)
        end


        parent_sub, dF = next_parent_dF(ppm, swaps, new_substitutions, delta_F, lineage_habit; parent = parent_sub)
        parent_fitness += dF
        # advance lineage


        fitness_log[gen + 1] = parent_fitness

        if !solved
            if parent_sub == inv_target
                solve_num = gen
                solved = true
            end
        end
    end


    return fitness_log, solve_num
end