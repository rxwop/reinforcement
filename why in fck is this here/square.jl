using RxCiphers
import Base.show



function vchecksquare(vector::Vector) ::Int
    L = length(vector)
    n = isqrt(L)
    
    if n ^ 2 != L
        e = ErrorException("Vector is not a square length")
        throw(e)
    end

    return n
end

mutable struct SquareKey
    charspace::NCharSpace{1}
    key::Vector{Int}
    size::Int

    displaykey::Matrix{String}

    tosquare::Dict{Int, Int}
    fromsquare::Dict{Int, Int}

    function SquareKey(key::Vector{Int}, charspace::NCharSpace{1}, merged::NTuple{N, Pair{Int, Int}} where N)
        # check that merged is good


        if !allunique(key)
            e = ArgumentError("All square key entries must be unique")
            throw(e)
        end

        n = vchecksquare(key)
        display = fill("", (n, n))

        tosquare = Dict{Int, Int}()
        fromsquare = Dict{Int, Int}()
        for squaretoken in 1:n^2
            exttoken = key[squaretoken]
            fromsquare[squaretoken] = exttoken
            tosquare[exttoken] = squaretoken

            display[squaretoken] = charspace.charmap[exttoken]
        end
        for (mergedtoken, mergedest) in merged
            squaretoken = tosquare[mergedest]
            tosquare[mergedtoken] = squaretoken

            display[squaretoken] *= "/$(charspace.charmap[mergedtoken])" 
        end


        return new(charspace, key, n, display, tosquare, fromsquare)
    end
end

# if merged not provided
# auto merge excluded tokens
function SquareKey(key::Vector{Int}, charspace::NCharSpace{1}) ::SquareKey
    # check it's not ridiculous to merge
    if charspace.size > 2 * length(key)
        e = ArgumentError("Too many tokens remain to be merged, choose a larger square")
        throw(e)
    end

    # you know exactly how long the vector should be, preallocate.
    tomerge = Vector{Int}()
    for i in 1:charspace.size
        if !(i in key)
            push!(tomerge, i)
        end
    end

    merged = Tuple(i => j for (i, j) in zip(tomerge, key))

    return SquareKey(key, charspace, merged)
end

SquareKey(txt::Txt, merged::Vararg{Pair{Int, Int}, N} where N) = SquareKey(checktokenised(txt), txt.charspace, merged)


function show(io::IO, S::SquareKey)
    show(io, S.displaykey)
end

function show(io::IO, ::MIME"text/plain", S::SquareKey)
    println(io, "$(S.size)-Square Key:")
    display(S.displaykey)
end




function mapsquare(vect::Vector{Int}, square::SquareKey)
    return [square.tosquare[i] for i in vect]
end

function invmapsquare(vect::Vector{Int}, square::SquareKey)
    return [square.fromsquare[i] for i in vect]
end