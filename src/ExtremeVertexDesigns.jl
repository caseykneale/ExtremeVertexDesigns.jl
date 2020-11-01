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
        Note: k is a ridge parameter to invert a matrix. This will be swapped out for an
        invertible starting design algorithm later...
    """
    function federov_D(Ndesign::Int, candidates::Matrix; k = 1e-6)
        Nsupport, V = size(candidates)
        @assert Ndesign < Nsupport "Design points must be < than support candidates."
        init_samples = sample(1:Nsupport, Ndesign; replace = false)
        Xdesign = candidates[init_samples,:]
        Xsupport = candidates[ setdiff(1:Nsupport, init_samples),:]
        Nsupport -= Ndesign 
        #ToDo: implement invertible starting design algorithm.
        tihkinov = Diagonal(ones(V) .* k)
        
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

end

