module QAOApools
    
    import PauliOperators
    import PauliOperators: FixedPhasePauli, Pauli, ScaledPauli
    import PauliOperators: ScaledPauliVector, PauliSum
    
    """
        qaoa_mixer(n::Int64)

    Returns the pool containing only the standard qaoa mixer.

    # Parameters
    - `n`: Number of qubits

    # Returns
    - `pool`: one element pool with qaoa mixer
    """
    function qaoa_mixer(n::Int64)
        pool = ScaledPauliVector{n}[]
        
        term = []
        print(term)
        for i in 1:n
            push!(term, ScaledPauli(Pauli(n; X=i)))
        end

        push!(pool, term)
        return pool
    end

    """
        qaoa_single_x(n::Int64)
    
    Returns the pool containing single-qubit Pauli Xs.

    # Parameters
    - `n`: Number of qubits
    
    # Returns
    - `pool`: pool containing Xs only
    """
    function qaoa_single_x(n::Int64)
        pool = ScaledPauliVector{n}[]
        
        for i in 1:n
            push!(pool, [ScaledPauli(Pauli(n; X=i))])
        end
        return pool
    end

    """
        qaoa_single_pool(n::Int64)
    
    Returns the pool containing single-qubit Pauli Xs and standard mixer.

    # Parameters
    - `n`: Number of qubits
    
    # Returns
    - `pool`: single-qubit pool
    """
    function qaoa_single_pool(n::Int64)
        return vcat(
            qaoa_single_x(n),
            qaoa_mixer(n)
        ) 
    end

    """
        qaoa_double_ops(n::Int64)
    
    Returns the pool containing two-qubit Paulis respecting bit-flip symmetry.

    # Parameters
    - `n`: Number of qubits
    
    # Returns
    - `pool`: pool containing symmetric two-qubit Paulis
    """
    function qaoa_double_ops(n::Int64)
        pool = ScaledPauliVector{n}[]

        for i in 1:(n-1)
            for j in (i+1):n
                push!(pool, [ScaledPauli(Pauli(n; X=[i,j]))])
                push!(pool, [ScaledPauli(Pauli(n; Y=[i,j]))])
                push!(pool, [ScaledPauli(Pauli(n; Y=i, Z=j))])
                push!(pool, [ScaledPauli(Pauli(n; Z=i, Y=j))])
            end
        end
        return pool
    end

    """
        qaoa_double_pool(n::Int64)
    
    Returns the pool containing symmetric single- and double-qubit Paulis and standard mixer.

    # Parameters
    - `n`: Number of qubits
    
    # Returns
    - `pool`: double-qubit pool
    """
    function qaoa_double_pool(n::Int64)
        return vcat(
            qaoa_single_x(n),
            qaoa_mixer(n),
            qaoa_double_ops(n)
        ) 
    end

end
