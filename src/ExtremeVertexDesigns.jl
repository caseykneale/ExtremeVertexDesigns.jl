module ExtremeVertexDesigns

    using Polyhedra, CDDLib, JuMP, Statistics, StatsBase, LinearAlgebra

    """
    extreme_values( model::JuMP.Model; barycenters = true)::Matrix{Float64}

        From an inequality constrained JuMP model, return a matrix containing candidate 
        design points.
        Note: Only barycenters of faces are currently supported. Soon, interior barycenters will also be
        added.
    """
    function extreme_values( model::JuMP.Model; barycenters = true)::Matrix{Float64}
        #ToDo: Invite other backends?
        poly = polyhedron(model, CDDLib.Library())
        removehredundancy!(poly)
        #ToDo: Faster memcpy then vcat from pointer?
        candidates = transpose(hcat(points(poly)...))
        if barycenters
            barycenters = [ mean( incidentpoints(poly, hidx) ) for hidx in eachindex(halfspaces(poly))]
            #ToDo: Faster memcpy then vcat from pointer?
            barycenters = transpose( hcat( barycenters...) )
            grand_center = mean( barycenters, dims = 1 )
            candidates = vcat( candidates, barycenters, grand_center )
        end
        return candidates
    end
    export extreme_values

    """
    federov_D(Ndesign::Int, candidates::Matrix; k = 1e-6)

        Performs Federov exchange algorithm to select a D-optimal design from a candidate matrix.
        Follows: Cook & Nachtsheim (1980).
        Note: λ is a ridge parameter to invert a matrix. This will be swapped out for an
        invertible starting design algorithm later...
    """
    function federov_D(Ndesign::Int, candidates::Matrix; λ = 1e-6)
        Nsupport, V = size(candidates)
        @assert Ndesign < Nsupport "Design points must be < than support candidates."
        init_samples = sample(1:Nsupport, Ndesign; replace = false)
        Xdesign = candidates[init_samples,:]
        Xsupport = candidates[ setdiff(1:Nsupport, init_samples),:]
        Nsupport -= Ndesign 
        #ToDo: implement invertible starting design algorithm.
        tihkinov = Diagonal(ones(V) .* λ)
        
        while(true) 
            XdtXd       = inv( Xdesign' * Xdesign + tihkinov )
            σ_design    = argmax( [ r' * XdtXd * r for r in eachrow( Xdesign ) ] )
            σ_support   = argmax( [ r' * XdtXd * r for r in eachrow( Xsupport ) ] )
            biggest_ds  = -Inf
            pair        = (0,0)
            for nd in 1:Ndesign, ns in 1:Nsupport
                σ_ds    = Xdesign[nd,:]' * XdtXd * Xsupport[ns,:]
                ds      = σ_support - ((σ_design * σ_support) - (σ_ds ^ 2) ) - σ_design
                if ds > biggest_ds
                    biggest_ds = ds
                    pair    = (nd, ns)
                end
            end
            if biggest_ds <= 0
                break
            end
            tmp = Xdesign[first(pair),:]
            Xdesign[first(pair),:] = Xsupport[last(pair),:]
            Xsupport[last(pair),:] = tmp
        end
        return (    score = det(inv( Xdesign' * Xdesign + tihkinov )), 
                    design = Xdesign )
    end
    export federov_D

    """
    K_exchange_D(k::Int, Ndesign::Int, candidates::Matrix; ridge = 1e-6)

        Performs the k-exchange algorithm to select a D-optimal design from a candidate matrix.
        Follows: Johnson & Nachtsheim (1983).
        Note: λ is a ridge parameter to invert a matrix. This will be swapped out for an
        invertible starting design algorithm later...
    """
    function k_exchange_D(k::Int, Ndesign::Int, candidates::Matrix; λ = 1e-6)
        Nsupport, V = size(candidates)
        @assert Ndesign < Nsupport "Design points must be < than support candidates."
        @assert k < Ndesign "Cannot exchange more rows(k) than are present in the design."
        init_samples = sample(1:Nsupport, Ndesign; replace = false)
        Xdesign = candidates[init_samples,:]
        Xsupport = candidates[ setdiff(1:Nsupport, init_samples),:]
        Nsupport -= Ndesign 
        #ToDo: implement invertible starting design algorithm.
        tihkinov = Diagonal(ones(V) .* λ)
        
        while(true) 
            XdtXd       = inv( Xdesign' * Xdesign + tihkinov )
            σ_design    = [ r' * XdtXd * r for r in eachrow( Xdesign ) ]
            biggest_ds = -Inf
            for nd in sortperm( σ_design )[1:k]
                biggest_ds = -Inf
                pair = (0,0)
                σ_k_design  = Xdesign[nd,:]' * XdtXd * Xdesign[nd,:]
                for ns in 1:Nsupport    
                    σ_support   = Xsupport[ns,:]' * XdtXd * Xsupport[ns,:]
                    σ_ds    = Xdesign[nd,:]' * XdtXd * Xsupport[ns,:]
                    ds      = σ_support - ((σ_k_design * σ_support) - (σ_ds ^ 2) ) - σ_k_design
                    if ds > biggest_ds
                        biggest_ds = ds
                        pair    = (nd, ns)
                    end
                end
                (biggest_ds <= 0) && break
                tmp = Xdesign[first(pair),:]
                Xdesign[first(pair),:] = Xsupport[last(pair),:]
                Xsupport[last(pair),:] = tmp    
            end
            (biggest_ds <= 0) && break
        end
        return (    score = det(inv( Xdesign' * Xdesign + tihkinov )), 
                    design = Xdesign )
    end
    export k_exchange_D

    """
        SimplexLatticeDesign( Components::Int, Spaces::Int )
    Returns an array of tuples (`Components` in length) which represent design points.
    """
    function SimplexLatticeDesign( Components::Int, Spaces::Int )
        @assert( ( Components > 0 ) && ( Spaces > 0 ), "Components and Spaces must be positive integers!" )
        DesignPoints = Inf
        try
            DesignPoints = Int( factorial( Components + Spaces - 1) / ( factorial(Spaces) * factorial(Components - 1) ))
        catch
            @warn("Cannot estimate number of design points. May execute for the lifetime of the universe.")
        end
        PossiblePoints = reverse(0:Spaces) ./ Spaces
        DesignSpace = repeat([PossiblePoints], Components)
        Possible = [ prod for prod in Base.product(DesignSpace...) if (sum(prod) == 1) ]

        if !isinf(DesignPoints)
            @assert( DesignPoints == length(Possible) )
        end
        return Possible
    end
    return SimplexLatticeDesign

    """
        SimplexCentroidDesign( Components::Int, Order::Union{ UnitRange{ Int }, Int} )
    Returns an array of tuples (`Components` in length) which represent design points.
    Please note: the Order can be an Int or a UnitRange.
    """
    function SimplexCentroidDesign( Components::Int, Order::Union{ UnitRange{ Int }, Int} )
        HighestOrder = Order[end]
        @assert( ( Components > 0 ) && ( HighestOrder > 0 ), "Components and the highest orders must be positive integers!" )
        @assert( HighestOrder <= Components, "The highest order of the design must be less than the number of components." )
        DesignPoints = Inf
        try
            DesignPoints = 2 ^ ( Components ) - 2 ^ ( Components - HighestOrder )
            if isa(Order, Int)
                DesignPoints = 2 ^ ( HighestOrder ) - 1
            end
        catch
            @warn("Cannot estimate number of design points. May execute for the lifetime of the universe.")
        end

        PossiblePoints = zeros( Components )
        DesignSpace = []
        for order in Order
            PossiblePoints[ 1:order ] .= 1.0 / order
            push!( DesignSpace, unique( collect( Combinatorics.permutations(PossiblePoints, Components) ) ) )
        end
        DesignSpace = ( isa(Order, Int) ) ? DesignSpace[1] : reduce(vcat, DesignSpace)
        if !isinf(DesignPoints)
            println( DesignPoints )
            println(length(DesignSpace))
            @assert( DesignPoints == length(DesignSpace) )
        end
        return DesignSpace
    end
    return SimplexCentroidDesign

end

