#= Set up various operator pools. =#
#= TO DO: do the extra functions in MyPauliOperators need to be added to the PauliOperators.jl package? =#
#= TO DO: tile operators by removing leading and trailing I's first. =#

module Pools

    import PauliOperators
    import PauliOperators: FixedPhasePauli, Pauli, ScaledPauli
    import PauliOperators: ScaledPauliVector, PauliSum

    import ..MyPauliOperators: otimes, ≈
    import PauliOperators: ⊗, ≈

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
        qubitexcitationpool(n_system::Int)
            
    The number of singles excitations = (n 2), and the doubles = 3*(n 4).
            
    # Parameters
    - `n_system`: Number of qubits in the system

    # Returns
    - `pool`: the qubit-excitation-based pool as defined in Communications Physics 4, 1 (2021).
    - `target_and_source`: Dict mapping each pool operator to the target and source orbitals involved in the excitation. 
    """               
    function qubitexcitationpool(n_system::Int)   
        pool = ScaledPauliVector{n_system}[]
        target_and_source = Dict{ScaledPauliVector{n_system}, Vector{Vector{Int64}}}()
                            
        for i in 1:n_system
            for j in i+1:n_system
                # singles excitations
                op = qubitexcitation(n_system, i, j)
                push!(pool, op)
                target_and_source[op] = [[i,j]]

                # doubles excitations
                for k in j+1:n_system
                    for l in k+1:n_system
                        target_pair = [i,j]; source_pair = [k,l]
                        new_op = qubitexcitation(n_system, target_pair[1], target_pair[2], source_pair[1], source_pair[2])
                        push!(pool, new_op)
                        target_and_source[new_op] = [target_pair,source_pair]

                        target_pair = [i,k]; source_pair = [j,l]
                        new_op = qubitexcitation(n_system, target_pair[1], target_pair[2], source_pair[1], source_pair[2])
                        push!(pool, new_op)
                        target_and_source[new_op] = [target_pair,source_pair]

                        target_pair = [j,k]; source_pair = [i,l]
                        new_op = qubitexcitation(n_system, target_pair[1], target_pair[2], source_pair[1], source_pair[2])
                        push!(pool, new_op)
                        target_and_source[new_op] = [target_pair,source_pair]
                    end
                end
            end
        end
        return pool, target_and_source                                            
    end

                                
    function qubitexcitation_complemented(n::Int, i::Int, k::Int)
        return (1/2) .* [
             ScaledPauli(Pauli(n; X=[i,k])),
             ScaledPauli(Pauli(n; Y=[i,k])),
        ]
    end
    function qubitexcitation_complemented(n::Int, i::Int, j::Int, k::Int, l::Int)
        return (1/8) .* [
             -ScaledPauli(Pauli(n; X=[i,j,k,l])),
             ScaledPauli(Pauli(n; X=[i,j], Y=[k,l])),
             -ScaledPauli(Pauli(n; X=[i,k], Y=[j,l])),
             -ScaledPauli(Pauli(n; X=[i,l], Y=[j,k])),
            -ScaledPauli(Pauli(n; X=[j,k], Y=[i,l])),
            -ScaledPauli(Pauli(n; X=[j,l], Y=[i,k])),
            ScaledPauli(Pauli(n; X=[k,l], Y=[i,j])),
            -ScaledPauli(Pauli(n; Y=[i,j,k,l])),
        ]
    end

    """
        qubitexcitationpool_complemented(n_system::Int)
                                
    Returns the complemented qubit excitation pool on n_system qubits, inspired from arXiv 2109.01318.

    # Parameters
    - `n_system`: Number of qubits in the system

    # Returns
    - `pool`: the complemented qubit-excitation-based pool.
    - `target_and_source`: Dict mapping each pool operator to the target and source orbitals involved in the excitation. 
    """                                
    function qubitexcitationpool_complemented(n_system::Int)
        pool = ScaledPauliVector{n_system}[]
        target_and_source = Dict{ScaledPauliVector{n_system}, Vector{Vector{Int64}}}()
                                            
        for i in 1:n_system
            for j in i+1:n_system
                # singles excitations
                op = qubitexcitation_complemented(n_system, i, j)
                push!(pool, op)
                target_and_source[op] = [[i,j]]

                # doubles excitations
                for k in j+1:n_system
                    for l in k+1:n_system
                        target_pair = [i,j]; source_pair = [k,l]
                        new_op = qubitexcitation_complemented(n_system, target_pair[1], target_pair[2], source_pair[1], source_pair[2])
                        push!(pool, new_op)
                        target_and_source[new_op] = [target_pair,source_pair]

                        target_pair = [i,k]; source_pair = [j,l]
                        new_op = qubitexcitation_complemented(n_system, target_pair[1], target_pair[2], source_pair[1], source_pair[2])
                        push!(pool, new_op)
                        target_and_source[new_op] = [target_pair,source_pair]

                        target_pair = [j,k]; source_pair = [i,l]
                        new_op = qubitexcitation_complemented(n_system, target_pair[1], target_pair[2], source_pair[1], source_pair[2])
                        push!(pool, new_op)
                        target_and_source[new_op] = [target_pair,source_pair]
                    end
                end
            end
        end
        return pool, target_and_source 
    end

    """
        qubitadaptpool(n_system::Int)
                                            
    Returns the qubit ADAPT pool on n_system qubits as defined in PRX QUANTUM 2, 020310 (2021). 
    It is generated by taking each qubit-excitation-based operator and breaking it into individual Pauli terms.

    # Parameters
    - `n_system`: Number of qubits in the system

    # Returns
    - `pool`: the qubit-ADAPT pool.                                                
    """                                     
    function qubitadaptpool(n_system::Int)
        pool = ScaledPauliVector{n_system}[]
                                                                
        for i in 1:n_system
            for j in i+1:n_system
                # singles operators
                push!(pool, [ScaledPauli(Pauli(n_system; X=[i], Y=[j]))])
                push!(pool, [ScaledPauli(Pauli(n_system; X=[j], Y=[i]))])

                # doubles operators
                for k in j+1:n_system
                    for l in k+1:n_system
                        push!(pool, [ScaledPauli(Pauli(n_system; X=[i,k,l], Y=[j]))])
                        push!(pool, [ScaledPauli(Pauli(n_system; X=[j,k,l], Y=[i]))])
                        push!(pool, [ScaledPauli(Pauli(n_system; X=[l], Y=[i,j,k]))])
                        push!(pool, [ScaledPauli(Pauli(n_system; X=[k], Y=[i,j,l]))])
                        push!(pool, [ScaledPauli(Pauli(n_system; X=[i,j,l], Y=[k]))])
                        push!(pool, [ScaledPauli(Pauli(n_system; X=[i,j,k], Y=[l]))])
                        push!(pool, [ScaledPauli(Pauli(n_system; X=[j], Y=[i,k,l]))])
                        push!(pool, [ScaledPauli(Pauli(n_system; X=[i], Y=[j,k,l]))])
                    end
                end
            end
        end
        return pool
    end

    """                
        fullpauli(n::Int)
            
    The pool of all (4^n) n-qubit Pauli operators.

    # Parameters
    - `n`: Number of qubits in the system

    # Returns
    - `pool`: the full pauli pool.
    """ 
    function fullpauli(n::Int)
        pool = ScaledPauliVector{n}[]
        for plist in Iterators.product(ntuple(i->["I","X","Y","Z"],n)...)
            pstr = join(plist)
            pauli = [ScaledPauli(Pauli(pstr))]
            push!(pool, pauli)
        end
        return pool
    end

    """
        one_local_pool(n::Int64, axes=["I","X","Y","Z"])
                                            
    Returns the one-local pool containing each one-local operator on n qubits. 

    # Parameters
    - `n`: Number of qubits in the system

    # Returns
    - `pool`: the one-local pool.                                                
    """                                                                          
    function one_local_pool(n::Int64, axes=["I","X","Y","Z"])
        pool = ScaledPauliVector(n)
        for i in 1:n
            "X" in axes && (push!(pool, ScaledPauli(Pauli(n; X=i))))
            "Y" in axes && (push!(pool, ScaledPauli(Pauli(n; Y=i))))
            "Z" in axes && (push!(pool, ScaledPauli(Pauli(n; Z=i))))
        end
        return pool
    end
         
    """
        two_local_pool(n::Int64, axes=["X","Y","Z"])
                                            
    Returns the two-local pool containing each two-local operator on n qubits. 

    # Parameters
    - `n`: Number of qubits in the system

    # Returns
    - `pool`: the two-local pool.                                                
    """                                                                                     
    function two_local_pool(n::Int64, axes=["X","Y","Z"])
        pool = ScaledPauliVector{n}[]
        for pair in Iterators.product(ntuple(i->1:n, 2)...)
            i,j = pair
            if i < j
                for pair2 in Iterators.product(ntuple(i->axes, 2)...)
                    a,b = pair2
                    if a == "I" || b == "I" 
                        continue
                    end
                    l = "I"^(i-1)*a*("I"^(j-i-1))*b*"I"^(n-j)
                    pauli = [ScaledPauli(Pauli(l))]
                    push!(pool, pauli)
                end
            end
        end
        return pool
    end
                                                                                                
    """
        oneandtwo_local_pool(n::Int64)
                                            
    Returns the union of the one-local and two-local pools on n qubits. 

    # Parameters
    - `n`: Number of qubits in the system

    # Returns
    - `pool`: union of one-local and two-local pools.                                   
    """                                                                                                
    function oneandtwo_local_pool(n::Int64)
        return vcat(
            one_local_pool(n),
            two_local_pool(n),
        )
    end

    """
        minimal_complete_pool(n::Int64)

    Return the minimal complete pool on `n` qubits corresponding to 
    the `V` pool in the qubit-ADAPT paper (PRX QUANTUM 2, 020310 (2021)).
    """
    function minimal_complete_pool(n::Int64)
        @assert n > 1 "This pool is defined for n>=2 qubits."
        res = Array{String,1}(["ZY","YI"])
        for i=3:n
            tempres = Array{String,1}()
            for pstr in res
                push!(tempres,"Z"*pstr)
            end
            push!(tempres,"Y"*("I"^(i-1)))
            push!(tempres,"IY"*("I"^(i-2)))
            res = copy(tempres)
        end
        pool = [[ScaledPauli(Pauli(pstr))] for pstr in res]
        return pool
    end

    #= TODO: Where should this function live? =#
    """
        tile_operators(L1::Int, L2::Int, chosen_operators::Vector{Vector{ScaledPauli{N}}}, PBCs)
                                                                                                                    
    Constructs the tiled operators for a system of `L2` qubits, given a set of operators
    defined for a smaller problem instance on `L1` qubits.

    # Parameters
    - `L1`: number of qubits for small problem instance 
    - `L2`: number of qubits for large problem instance
    - `chosen_operators`: list of operators for small problem instance
    - `PBCs`: periodic boundary conditions

    # Returns
    - `tiled_ops`: tiled operators as a Vector{Vector{ScaledPauli}}
    """                                                                                                                    
    function tile_operators(L1::Int, L2::Int, chosen_operators::Vector{Vector{ScaledPauli{N}}}, PBCs) where {N}
        @assert L2 >= L1 "L2 must be greater than L1."
        if L2 == L1
            return chosen_operators
        end
        l = length(chosen_operators)
        tiled_ops = ScaledPauliVector{L2}[]
        permutations = PBCs ? (L2) : (L2-L1+1)
        for i=1:l
            spv = chosen_operators[i]
            for j = 1:permutations
                tiled_spv = ScaledPauli{L2}[]
                id_left = Pauli("I"^(j-1))
                id_right = Pauli("I"^(L2-L1+1-j))
                for sp in spv
                    tiled_sp = (j-1)==0 ? sp : otimes(id_left,sp)
                    tiled_sp = (L2-L1+1-j)==0 ? tiled_sp : otimes(tiled_sp, id_right)
                    push!(tiled_spv, tiled_sp)
                end
                duplicate = false
                for op in tiled_ops
                    if tiled_spv ≈ op
                        duplicate = true
                        break
                    end
                end
                if !duplicate
                    push!(tiled_ops, tiled_spv)
                end            
            end
        end
        return tiled_ops
    end
end
                                                                                                                                 