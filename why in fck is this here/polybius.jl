include("square.jl")


mutable struct Polybius <: AbstractCipher
    key::SquareKey
    indices::Vector{String}
    outcharspace::NCharSpace{2}

    inverted::Bool

    function Polybius(key::SquareKey, indices::Vector{String} = ["1", "2", "3", "4", "5"])
        W = CharSpace(indices)

        return new(key, indices, W ^ 2, false)
    end
end

invert!(P::Polybius) = RxCiphers.switch_invert_tag!(P)

function apply!(P::Polybius, txt::Txt) ::Txt
    if !txt.is_tokenised
        throw(TokeniseError)
    end

    if P.inverted
        # take text
        # nchar
        nchar!(txt, 2)

        if txt.charspace != P.outcharspace
            println("WARNING: Txt character space does not match Polybius output space")
        end

        # reassign to key.inspace
        txt.charspace = key.charspace

        # map from square
        txt.tokenised = invmapsquare(txt.tokenised, P.key)

        # bye bye
    else
        if txt.charspace != P.key.charspace
            println("WARNING: Txt character space does not match Polybius input space")
        end

        # take text
        # map to square
        txt.tokenised = mapsquare(txt.tokenised, P.key)

        # reassign to NCharSpace{2}
        txt.charspace = P.outcharspace

        # reduce
        reduce!(txt)

        # bye bye
    end

    return txt
end