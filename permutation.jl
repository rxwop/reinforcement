include("charspace.jl")
include("cipher.jl")
include("array functions.jl")
import Base.show
using Combinatorics


mutable struct Permutation <: AbstractCipher
    n::Int
    permutation::Vector{Int}

    remove_nulls::Bool
    inverted::Bool

    # UNSAFE DO NOT USE
    function Permutation(n::Int, permutation::Vector{Int}, remove_nulls::Bool, inverted::Bool)
        new(n, checkperm(permutation), remove_nulls, inverted)
    end
    # UNSAFE DO NOT USE


end
Permutation(permutation::Vector{Int}, remove_nulls::Bool = false) = Permutation(length(permutation), permutation, remove_nulls, false)
Permutation(n::Int, remove_nulls::Bool = false) = Permutation(n, collect(1:n), remove_nulls, false)
invPermutation(permutation::Vector{Int}, remove_nulls::Bool = false) = Permutation(length(permutation), permutation, remove_nulls, true)
invPermutation(n::Int, remove_nulls::Bool = false) = Permutation(n, collect(1:n), remove_nulls, true)

function invert!(P::Permutation)
    switch_invert_tag!(P)
    P.permutation = invperm(P.permutation)
    return P
end


function switch(P::Permutation, posa::Int, posb::Int) ::Permutation
    return Permutation(P.n, switch(P.permutation, posa, posb), P.remove_nulls, P.inverted)
end

function switch!(P::Permutation, posa::Int, posb::Int) ::Permutation
    switch!(P.permutation, posa, posb)
    return P
end

function splice!(P::Permutation, ind::Int, start::Int, finish::Int) ::Permutation
    splice!(P.permutation, ind, start, finish)
    return P
end

function splice(P::Permutation, ind::Int, start::Int, finish::Int) ::Permutation
    return Permutation(P.n, splice(P.permutation, ind, start, finish), P.remove_nulls, P.inverted)
end


function show(io::IO, T::Permutation)
    show(io, T.permutation)
end

function show(io::IO, ::MIME"text/plain", T::Permutation)
    if T.inverted
        inverse_text = " (padding inverted)"
    else
        inverse_text = ""
    end

    println(io, "$(T.n)-column Permutation", inverse_text, ":")
    show(io, T.permutation)
end




function reinsert_nulls(vect::Vector{Int}, T::Permutation) ::Vector{Int}
    L = length(vect)
    overhang = L % T.n

    if overhang == 0
        return vect
    end

    vector = copy(vect)

    null_ends = T.permutation .> overhang

    num_full_rows = floor(Int, L / T.n)

    for (i,j) in enumerate(null_ends)
        if j
            insert!(vector, num_full_rows * T.n + i, NULL_TOKEN)
        end
    end

    return vector
end



function unshape(vect::Vector{Int}, T::Permutation) ::Matrix{Int}
    new_tokens = reinsert_nulls(vect, T)
    new_tokens = reshape(new_tokens, (T.n, :))

    return new_tokens
end



function apply(T::Permutation, vect::Vector{Int}; safety_checks::Txt) ::Vector{Int}
    if T.inverted # inverse application
        new_tokens = unshape(vect, T)
        new_tokens = vec(new_tokens[T.permutation, :])

        if T.remove_nulls
            return filter!(!=(NULL_TOKEN), new_tokens)
        else
            return new_tokens
        end

    else # regular application
        new_tokens = safe_reshape_2D(vect, T.n, NULL_TOKEN)
        new_tokens = new_tokens[T.permutation, :]

        new_tokens = vec(new_tokens)

        if T.remove_nulls
            return filter!(!=(NULL_TOKEN), new_tokens)
        else
            return new_tokens
        end

    end

end