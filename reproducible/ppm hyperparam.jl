include("ppm keysearch.jl")

using IterTools
iter_M = 10:10:30
iter_spawns = 2:4:14
iter = IterTools.product(iter_M, iter_spawns)

m_size = size(iter)

percentages = Matrix{Float64}(undef, m_size)
key_num = Matrix{Float64}(undef, m_size)

for (i, (M, spawns)) in enumerate(iter)
    (k, p) = avg_keys_searched(75, 75, spawns, float(M))
    percentages[i] = p
    key_num[i] = k
end


########## EXPORT DATA ###########


using Plots.PlotMeasures, Plots

p1 = heatmap(iter_spawns, iter_M, key_num,
    color = cgrad(:acton),
    aspect_ratio = :equal,
    colorbar = false,
    xticks = iter_spawns,
    yticks = iter_M,
)

for (i, (M, spawns)) in enumerate(iter)
    annotate!(p1,
        spawns,
        M,
        text(round(key_num[i]; digits = 2); pointsize = 8),
    )
end

ylabel!(p1, "Update Rate")
xlabel!(p1, "Spawn Number")
title!(p1, "Keys Searched")

p2 = heatmap(iter_spawns, iter_M, percentages,
    color = cgrad(:acton),
    aspect_ratio = :equal,
    colorbar = false,
    xticks = iter_spawns,
    yticks = iter_M,
)

for (i, (M, spawns)) in enumerate(iter)
    annotate!(p2,
        spawns,
        M,
        text(round(percentages[i]; digits = 2); pointsize = 8),
    )
end

ylabel!(p2, "Update Rate")
xlabel!(p2, "Spawn Number")
title!(p2, "Percentage Solved")

p = plot(p1, p2,
margins = 20px,
thickness_scaling = 0.7
)
