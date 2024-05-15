import LinearAlgebra: I

import PauliOperators
import PauliOperators: FixedPhasePauli, Pauli, ScaledPauli
import PauliOperators: ScaledPauliVector, PauliSum
import PauliOperators: clip!, jordan_wigner

"""
    xyz_model(L::Int, Jx::Float, Jy::Float, Jz::Float, PBCs::Bool)

An XYZ Heisenberg Hamiltonian.

# Parameters
- `L`: system size.
- `Jx`: coupling along X.
- `Jy`: coupling along Y.
- `Jz`: coupling along Z.
- `PBCs`: Periodic Boundary Conditions

# Returns
- `PauliOperators.PauliSum`: the Hamiltonian

"""
function xyz_model(L,Jx,Jy,Jz,PBCs)
    H = PauliSum(L)

    for i=1:L-1
        term = PauliSum(L)
        term += Jx*Pauli(L; X=[i,i+1])
        term += Jy*Pauli(L; Y=[i,i+1])        
        term += Jz*Pauli(L; Z=[i,i+1])        
        clip!(term)
        sum!(H,term)
    end
    if PBCs
        term = PauliSum(L)
        term += Jx*Pauli(L; X=[L,1])
        term += Jy*Pauli(L; Y=[L,1])        
        term += Jz*Pauli(L; Z=[L,1])        
        clip!(term)
        sum!(H,term)
    end

    return H
end

"""
    hubbard_jw(graph::Array{T,2}, U, t)

A Hubbard Hamiltonian in the Jordan-Wigner basis.

Copied shamelessly from Diksha's ACSE repository.

# Parameters
- `graph`: an adjacency matrix identifying couplings. Must be symmetric.
- `U`: Coulomb interaction for all sites
- `t`: hopping energy for all couplings

# Returns
- `PauliOperators.PauliSum`: the Hamiltonian

"""
function hubbard_hamiltonian(graph::Matrix{T}, U, t) where T
    Ni, Nj = size(graph)
    Ni == Nj || throw(DimensionMismatch)
    Norb = Ni
    N = 2*Norb
    H = PauliSum(N)
    for i in 1:Norb
        ia = 2*i-1
        ib = 2*i
        for j in i+1:Norb
            ja = 2*j-1
            jb = 2*j
            abs(graph[i,j]) > 1e-16 || continue

            tij = jordan_wigner(ia, N)*jordan_wigner(ja,N)'
            tij += jordan_wigner(ib, N)*jordan_wigner(jb,N)'
            tij += adjoint(tij)
            clip!(tij)
            sum!(H,t*tij)
        end

        ni = jordan_wigner(ia, N)*jordan_wigner(ia,N)'*jordan_wigner(ib, N)*jordan_wigner(ib,N)'
        sum!(H,U*ni)
    end
    return H
end

"""
    hubbard_hamiltonian(L::Int, U, t; pbc=false)

Convenience constructor for a 1D nearest-neighbor Hubbard model with L sites.

"""
function hubbard_hamiltonian(L::Int, U, t; pbc=false)
    A = Matrix{Bool}(I, L, L)
    pbc && (A[L,1] = A[1,L] = true)
    return hubbard_hamiltonian(A, U, t)
end
