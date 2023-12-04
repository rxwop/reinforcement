using RxCiphers
include("square.jl")

const fixed_playfairmap = [0, 28, 29, 30, 26, 136, 127, 128, 129, 130, 141, 252, 253, 254, 255, 146, 377, 378, 379, 380, 126, 502, 503, 504, 505, 52, 0, 54, 55, 51, 151, 162, 153, 154, 155, 276, 167, 278, 279, 280, 401, 172, 403, 
404, 405, 526, 152, 528, 529, 530, 77, 78, 0, 80, 76, 176, 177, 188, 179, 180, 301, 302, 193, 304, 305, 426, 427, 198, 429, 430, 551, 552, 178, 554, 555, 102, 103, 104, 0, 101, 201, 202, 203, 214, 205, 326, 327, 328, 219, 330, 451, 452, 453, 224, 455, 576, 577, 578, 204, 
580, 2, 3, 4, 5, 0, 226, 227, 228, 229, 240, 351, 352, 353, 354, 245, 476, 477, 478, 479, 250, 601, 602, 603, 604, 230, 256, 7, 8, 9, 10, 0, 158, 159, 160, 156, 266, 257, 258, 259, 260, 271, 382, 383, 384, 385, 251, 507, 508, 509, 510, 31, 282, 33, 34, 35, 182, 0, 184, 185, 181, 281, 292, 283, 284, 285, 406, 297, 408, 409, 410, 531, 277, 
533, 534, 535, 56, 57, 308, 59, 60, 207, 208, 0, 210, 206, 306, 307, 318, 309, 310, 431, 432, 323, 434, 435, 556, 557, 303, 559, 560, 81, 82, 83, 334, 85, 232, 233, 234, 0, 231, 331, 332, 333, 344, 335, 456, 457, 458, 349, 460, 581, 582, 583, 329, 585, 106, 107, 108, 109, 360, 132, 133, 134, 135, 0, 356, 357, 358, 359, 370, 481, 482, 483, 484, 375, 606, 607, 608, 609, 355, 381, 12, 13, 14, 15, 386, 137, 138, 139, 140, 0, 288, 289, 290, 286, 396, 387, 388, 389, 390, 376, 512, 513, 514, 515, 36, 407, 38, 39, 40, 161, 412, 163, 164, 165, 312, 0, 314, 315, 311, 411, 422, 413, 414, 415, 536, 402, 538, 539, 540, 61, 62, 433, 64, 65, 186, 187, 438, 189, 190, 337, 338, 0, 340, 336, 436, 437, 448, 439, 440, 561, 562, 428, 564, 565, 86, 87, 88, 459, 90, 211, 212, 213, 464, 215, 362, 363, 364, 0, 361, 461, 462, 463, 474, 465, 586, 587, 588, 454, 590, 111, 112, 113, 114, 485, 236, 237, 238, 239, 490, 262, 263, 264, 265, 0, 486, 487, 488, 489, 500, 611, 612, 613, 614, 480, 506, 17, 18, 19, 20, 511, 142, 143, 144, 145, 516, 267, 268, 269, 270, 0, 418, 419, 420, 416, 501, 517, 518, 519, 520, 41, 532, 43, 44, 45, 166, 537, 168, 169, 170, 291, 542, 293, 294, 295, 442, 0, 444, 445, 441, 541, 527, 543, 544, 545, 66, 67, 558, 69, 70, 191, 192, 563, 194, 195, 316, 317, 568, 319, 320, 467, 468, 0, 470, 466, 566, 567, 553, 569, 570, 91, 92, 93, 584, 95, 216, 217, 218, 589, 220, 341, 342, 343, 594, 345, 492, 493, 494, 0, 491, 591, 592, 593, 579, 595, 116, 117, 118, 119, 610, 241, 242, 243, 244, 
615, 366, 367, 368, 369, 620, 392, 393, 394, 395, 0, 616, 617, 618, 
619, 605, 6, 22, 23, 24, 25, 11, 147, 148, 149, 150, 16, 272, 273, 274, 275, 21, 397, 398, 399, 400, 0, 548, 549, 550, 546, 46, 32, 48, 
49, 50, 171, 37, 173, 174, 175, 296, 42, 298, 299, 300, 421, 47, 423, 424, 425, 572, 0, 574, 575, 571, 71, 72, 58, 74, 75, 196, 197, 63, 199, 200, 321, 322, 68, 324, 325, 446, 447, 73, 449, 450, 597, 598, 0, 600, 596, 96, 97, 98, 84, 100, 221, 222, 223, 89, 225, 346, 347, 348, 94, 350, 471, 472, 473, 99, 475, 622, 623, 624, 0, 621, 121, 122, 123, 124, 110, 246, 247, 248, 249, 115, 371, 372, 373, 374, 120, 496, 497, 498, 499, 125, 522, 523, 524, 525, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

const AtoY = CharSpace(collect("ABCDEFGHIJKLMNOPQRSTUVWXY"))^2

const BiAlph = Alphabet ^ 2

mutable struct Playfair <: AbstractCipher
    key::SquareKey

    splittoken::Int # confusingly with Alphabet tokenisation
    splitmap::Dict{Int, Tuple{Int, Int}}
    forbidden::Int

    fromsquaremap::Vector{Int}
    outcharspace::NCharSpace{2}

    inverted::Bool

    function Playfair(key::SquareKey, spl::Int)
        if length(key) != 25
            error("Input key must be of length 25 (omit tokens to be merged)")
        end

        W_base = CharSpace(key)
        W = W_base ^ 2
        fromsquaremap = Vector{Int}(undef, 25)
        for (tok, str) in enumerate(W_base.charmap)
            p = findfirst(==(str), Alphabet.charmap)
            fromsquaremap[tok] = p
        end

        forbiddenbigram = BiAlph.reducemap[spl, spl]

        splits = Dict{Int, Tuple{Int, Int}}()
        for i in 1:25
            isobigram = BiAlph.reducemap[i, i]
            splitbigram = (BiAlph.reducemap[i, spl], BiAlph.reducemap[spl, i])
            
            splits[isobigram] = splitbigram
        end

        return new(key, spl, splits, forbiddenbigram, fromsquaremap, W, false)
    end
end


function replace_insert(vector::Vector{Int}, dict::Dict)
    out = Vector{Int}()
    for i in vector
        append!(out, get(dict, i, i))
    end

    return out
end


function apply(P::Playfair, vector::Vector{Int}; safety_checks::Txt, playfairmap::Vector{Int} = fixed_playfairmap)
    # Split isobigrams
    # Apply fixed map

    if P.inverted
    else
        if P.forbidden in vector
            e = ErrorException("$(P.outcharspace.charmap[P.forbidden]) bigram cannot be split by $(P.outcharspace.charmap[P.splittoken])")
            throw(e)
        end

        new_tokens = replace_insert(vector, P.splitmap) # process in BiAlph
        new_tokens = [BiAlph.reducemap[i] for i in new_tokens] # return to monograms
        # reassign characters to square space
        # return to bigrams
        new_tokens = [playfairmap[i] for i in new_tokens] # execute map
        # return to monograms
        new_tokens = [P.fromsquaremap[i] for i in new_tokens] # return to Alphabet
    end


    return new_tokens
end







# cart25 = CartesianIndices(Array{Int}(undef, (5,5)))
# linear25 = LinearIndices(Array{Int}(undef, (5,5)))

# fixed_playfairmap = zeros(Int, 675)
# for (bigram, str) in enumerate(AtoY.charmap)
#     (tokenA, tokenB) = AtoY.reducemap[bigram]
#     if tokenA == tokenB
#         continue
#     end

#     A = cart25[tokenA]
#     B = cart25[tokenB]

#     if A[1] == B[1]

#         # increment
#         A += CartesianIndex(0, 1)
#         B += CartesianIndex(0, 1)

#         # wrap around
#         if B[2] == 6
#             B = CartesianIndex(B[1], 1)
#         elseif A[2] == 6
#             A = CartesianIndex(A[1], 1)
#         end

#     elseif A[2] == B[2]

#         # increment
#         A += CartesianIndex(1, 0)
#         B += CartesianIndex(1, 0)

#         # wrap around
#         if B[1] == 6
#             B = CartesianIndex(1, B[2])
#         elseif A[1] == 6
#             A = CartesianIndex(1, A[2])
#         end

#     else
#         # switch corners
#         a2, b2 = B[2], A[2]
#         A = CartesianIndex(A[1], a2)
#         B = CartesianIndex(B[1], b2)
#     end

#     tokenA = linear25[A]
#     tokenB = linear25[B]

#     newbigram = AtoY.reducemap[tokenA, tokenB]

#     fixed_playfairmap[bigram] = newbigram
# end
