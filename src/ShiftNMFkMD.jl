module ShiftNMFk

    using JSON
    using Clustering
    using Distances
    using Stats
    using Gadfly
    using Compose
    using NLopt

    export shiftNMFk
    export Plot
    export ResultsForNumSources
    export FindLocations
    export ParseLoc
    export AIC_final

    export Parallel_Tri
    export PTri
    export Triangulate

    include("ShiftNMF2.jl");
    include("Parallel_ShiftNMF.jl");   
    include("ParseTrials.jl");
    include("Plot.jl");
    include("ShiftNMFk_CosClust.jl");
    include("ShiftNMFk_using_functionsCos.jl");

    include("TriangulatePos.jl");
    include("Parallel_Triangle.jl");
    include("LocCluster.jl");

    include("ParseLoc.jl");
    include("AIC_final.jl");

end # module