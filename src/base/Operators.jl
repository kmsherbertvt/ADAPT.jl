"""
    Operators

A suite of common operators, especially useful for constructing operator pools.

TODO:
    I haven't decided yet whether observables should live here or not.
    If they do, I'll want to standardize the interface somehow.
    In particular, the interface with pyscf for molecules is rather hazy.
    I think we need a separate package which is a Julia wrapper for openfermion.
    Then observables will generally be input as qubit operators from that package,
        or perhaps we have a simple method that converts qubit operators to PauliSums,
        so we have better control over the arithmetic being performed.
    In any case, though I may evict them someday,
        standard lattice systems like Hubbard and Heisenberg, not requiring openfermion,
        may inhabit this module for the time being.

"""
module Operators
    import LinearAlgebra: I

    import PauliOperators
    import PauliOperators: FixedPhasePauli, Pauli, ScaledPauli
    import PauliOperators: ScaledPauliVector, PauliSum
    import PauliOperators: clip!, jordan_wigner


    """
        qubitexcitation(n::Int, i::Int, k::Int)
        qubitexcitation(n::Int, i::Int, j::Int, k::Int, l::Int)

    Qubit excitation operators as defined in Yordanov et al. 2021.

    Note that Yordanov's unitaries are defined as `exp(iθG)` rather than `exp(-iθG)`,
        so variational parameters will be off by a sign.

    # Parameters
    - `n`: total number of qubits
    - `i,j,k,l`: qubit indices as defined in Yordanov's paper.

    # Returns
    - `PauliOperators.ScaledPauliVector`: the qubit excitation operator

        Note that all Pauli terms in any single qubit excitation operator commute,
            so the `ScaledPauliVector` representation is "safe".

    """
    function qubitexcitation(n::Int, i::Int, k::Int)
        return (1/2) .* [
             ScaledPauli(Pauli(n; X=[i], Y=[k])),
            -ScaledPauli(Pauli(n; X=[k], Y=[i])),
        ]
    end

    # TODO: Karunya says this is not the only 2-body QEB operator.
    function qubitexcitation(n::Int, i::Int, j::Int, k::Int, l::Int)
        return (1/8) .* [
             ScaledPauli(Pauli(n; X=[i,k,l], Y=[j])),
             ScaledPauli(Pauli(n; X=[j,k,l], Y=[i])),
             ScaledPauli(Pauli(n; X=[l], Y=[i,j,k])),
             ScaledPauli(Pauli(n; X=[k], Y=[i,j,l])),
            -ScaledPauli(Pauli(n; X=[i,j,l], Y=[k])),
            -ScaledPauli(Pauli(n; X=[i,j,k], Y=[l])),
            -ScaledPauli(Pauli(n; X=[j], Y=[i,k,l])),
            -ScaledPauli(Pauli(n; X=[i], Y=[j,k,l])),
        ]
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



end