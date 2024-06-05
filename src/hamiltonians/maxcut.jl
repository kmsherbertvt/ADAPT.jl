#= All the methods here are slightly modified from Sam's AdaptBarren code. =#
module MaxCut
    using Random
    using Erdos
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

    """
        max_cut_hamiltonian(n::Int, edges::Vector{Tuple{Int, Int, T}}) where T<:Real

    Return the max cut Hamiltonian acting on `n` qubits which corresponds to the edge
    set `edges`.

    # Examples

    Elements of the edge set look like `(1,2,5.0)` corresponding to an edge between
    vertices `1` and `2` with edge weight `5.0`.
    """
    function max_cut_hamiltonian(n::Int, edges::Vector{Tuple{Int, Int, T}}) where T<:Real
        H = ScaledPauliVector{n}() 
        for (i,j,w)=edges
            term = ScaledPauliVector{n}() 
            term = (-w/2.0)*ScaledPauli(Pauli(n)) 
            push!(H, term)

            term = ScaledPauliVector{n}()
            term = (-w/2.0)*ScaledPauli(Pauli(n; Z=[i,j]))
            push!(H, term)
        end
        return H
    end

    function get_random_unweighted_graph_edges(n::Int, k::Int; rng = _DEFAULT_RNG)
        seed = abs(rand(rng, Int64) % 10^7)
        g = random_regular_graph(n, k, Network, seed=seed)
        return [(i,j,1.0) for (i,j) in Erdos.edges(g)]
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
        return max_cut_hamiltonian(n, v)
    end
end