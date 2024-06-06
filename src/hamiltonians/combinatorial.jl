using Graphs
using Random

import LinearAlgebra: I

import PauliOperators
import PauliOperators: FixedPhasePauli, Pauli, ScaledPauli
import PauliOperators: ScaledPauliVector, PauliSum
import PauliOperators: clip!, jordan_wigner

_DEFAULT_RNG = MersenneTwister(1234);

"""
    get_unweighted_maxcut(g::SimpleGraph)

Take a graph object and extract edges for MaxCut.

# Parameters
- `g`: graph instance.

# Returns
- `edge_list`: list of edges and weights equal to one.

"""
function get_unweighted_maxcut(g::SimpleGraph)
    edge_indices = edges(g)
    edge_list = [(src(e), dst(e), 1.0) for e in edge_indices]

    #for e in edge_indices
    #    push!(edge_list, (src(e), dst(e), 1.0))
    #end
    return edge_list
end


"""
    get_weighted_maxcut(g::SimpleGraph, rng = _DEFAULT_RNG)

Take a graph object and extract edges and assign edge weights.

# Parameters
- `g`: graph instance.
- `rng`: random number generator to generate weights.

# Returns
- `edge_list`: list of edges and weights.

"""
function get_weighted_maxcut(g::SimpleGraph, rng = _DEFAULT_RNG)
    edge_indices = edges(g)
    edge_list = [(src(e), dst(e), rand(rng, Float64)) for e in edge_indices]
    #for e in edge_indices
    #    push!(edge_list, (src(e), dst(e), rand(rng, Float64)))
    #end
    return edge_list
end


"""
    maxcut_hamiltonian(V::Int, Edges::Vector{Tuple{Int,Int,T}}) where T<:Real

A MaxCut Hamiltonian defined on a graph containing only Pauli ZZ terms.

# Parameters
- `V`: number of vertices.
- `Edges`: list of edges.

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

