import Graphs
import Random

import LinearAlgebra: I

import PauliOperators
import PauliOperators: FixedPhasePauli, Pauli, ScaledPauli
import PauliOperators: ScaledPauliVector, PauliSum
import PauliOperators: clip!, jordan_wigner

_DEFAULT_RNG = Random.MersenneTwister(1234);

"""
    get_unweighted_maxcut(g::Graphs.SimpleGraph)

Take a graph object and extract edges for MaxCut.

# Parameters
- `g`: graph instance.

# Returns
- `edge_list`: list of edges and weights equal to one.

"""
function get_unweighted_maxcut(g::Graphs.SimpleGraph)
    edge_indices = Graphs.edges(g)
    edge_list = [(Graphs.src(e), Graphs.dst(e), 1.0) for e in edge_indices]
    return edge_list
end


"""
    get_weighted_maxcut(g::Graphs.SimpleGraph, rng = _DEFAULT_RNG)

Take a graph object and extract edges and assign edge weights.

# Parameters
- `g`: graph instance.
- `rng`: random number generator to generate weights.

# Returns
- `edge_list`: list of edges and weights.

"""
function get_weighted_maxcut(g::Graphs.SimpleGraph, rng = _DEFAULT_RNG)
    edge_indices = Graphs.edges(g)
    edge_list = [(Graphs.src(e), Graphs.dst(e), rand(rng, Float64)) for e in edge_indices]
    return edge_list
end


"""
    maxcut_hamiltonian(V::Int, Edges::Vector{Tuple{Int,Int,T}}) where T<:Real

A MaxCut Hamiltonian defined on a graph containing only Pauli ZZ terms.

# Parameters
- `V`: number of vertices.
- `Edges`: list of edges, in the form of (first index, second index, weight).

# Returns
- `H`: MaxCut Hamiltonian

"""
function maxcut_hamiltonian(V::Int, Edges::Vector{Tuple{Int,Int,T}}) where T<:Real
    H = ScaledPauliVector{V}()

    for (i,j,w)=Edges
        term = (-w/2.0)*ScaledPauli(Pauli(V; Z=[i,j]))       
        push!(H,term)

        term = (-w/2.0)*ScaledPauli(Pauli(V)) 
        push!(H,term)
    end

    return H
end

