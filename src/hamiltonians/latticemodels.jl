module LatticeHamiltonians
    import PauliOperators: Pauli, PauliSum, clip!

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
end