import LinearAlgebra: I

import PauliOperators
import PauliOperators: FixedPhasePauli, Pauli, ScaledPauli
import PauliOperators: ScaledPauliVector, PauliSum
import PauliOperators: clip!, jordan_wigner

"""
    maxcut_hamiltonian(V::Int, Edges::Vector{Tuple{Int,Int,T}}) where T<Real

A MaxCut Hamiltonian defined on a graph containing only Pauli ZZ terms.

# Parameters
- `V`: number of vertices.
- `Edges`: list of edges.

# Returns
- `PauliOperators.PauliSum`: the Hamiltonian

"""
function maxcut_hamiltonian(V,Elist)
    H = PauliSum(V)

    for (i,j,w)=Edges
        term = PauliSum(V)
        term += (-w/2.0)*Pauli(V; Z=[i,j])       
        clip!(term)
        sum!(H,term)

        term = PauliSum(V)
        term += (-w/2.0)*Pauli(n)
        clip!(term)
        sum!(H,term)
    end

    return H
end

