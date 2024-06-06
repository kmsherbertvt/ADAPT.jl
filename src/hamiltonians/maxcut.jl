#= All the methods here are slightly modified from Sam's AdaptBarren code. =#
module MaxCut
    import Random: MersenneTwister
    import ..get_unweighted_maxcut  # function defined in combinatorial.jl
    import ..maxcut_hamiltonian  # function defined in combinatorial.jl
    import Graphs.SimpleGraphs.random_regular_graph
    import PauliOperators: FixedPhasePauli, Pauli, ScaledPauli, ScaledPauliVector, PauliSum
    import PauliOperators: clip!

    _DEFAULT_RNG = MersenneTwister(1234);

    function _id_x_str(n::Int, k::Int)
        s = repeat(["I"], n)
        s[k] = "X"
        return join(s)
    end

    function qaoa_mixer(n::Int)
        mixer = [1.0*Pauli(_id_x_str(n, k)) for k in range(1,n)] # ScaledPauliVector type
        return mixer
    end

    function get_random_unweighted_graph_edges(n::Int, k::Int; rng = _DEFAULT_RNG)
        g =  random_regular_graph(n,k)
        edge_list = get_unweighted_maxcut(g)
        return [(i,j,1.0) for (i,j) in edge_list]
    end

    function randomize_edge_weights!(v::Vector{Tuple{Int, Int, T}}; rng = _DEFAULT_RNG) where T<:Number
        for i=1:length(v)
            a, b, _ = v[i]
            v[i] = (a, b, rand(rng, Float64))
        end
    end


    """
        random_regular_max_cut_hamiltonian(n::Int, k::Int; rng = _DEFAULT_RNG, weighted = true)

    Return a random Hamiltonian for a max cut problem on `n` qubits.

    The corresponding graph is degree `k`. If an RNG is provided, this will be used to sample
    the graph and edge weights. If `weighted` is true, the edge weights will be randomly sampled
    from the uniform distribution `U(0,1)`.
    """
    function random_regular_max_cut_hamiltonian(n::Int, k::Int; rng = _DEFAULT_RNG, weighted = true)
        v = get_random_unweighted_graph_edges(n, k; rng=rng)
        if weighted
            randomize_edge_weights!(v; rng=rng)
        end
        return maxcut_hamiltonian(n, v)
    end
end