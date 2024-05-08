#= Set up various operator pools. =#

module OperatorPools
    import ..ADAPT
    import PauliOperators: Pauli, PauliSum, ScaledPauli, ScaledPauliVector, clip!, otimes, ⊗, get_phase, ≈
    import Combinatorics: combinations

    function fullpauli(n::Int)
        pool = ScaledPauliVector{n}[]
        for plist in Iterators.product(ntuple(i->["I","X","Y","Z"],n)...)
            pstr = join(plist)
            pauli = [ScaledPauli(Pauli(pstr))]
            push!(pool, pauli)
        end
    #    pool = pool[2:end] # skip the first entry, which is just "III...."  
        return pool
    end

    function qubitexcitationpool(n_system::Int)
        """
        The number of singles excitations = (n 2), and the doubles = 3*(n 4).
                
        qubitexcitationpool(n_system::Int)

        # Parameters
        - `n_system`: Number of qubits in the system
                
        # Returns
        - `pool`: the qubit-excitation-based pool as defined in Communications Physics 4, 1 (2021).
        - `target_and_source`: Dict mapping each pool operator to the target and source orbitals involved in the excitation. 
        """      
        pool = ScaledPauliVector{n_system}[]
        target_and_source = Dict{ScaledPauliVector{n_system}, Vector{Vector{Int64}}}()
        # single excitations
        for i in 1:n_system
            for j in i+1:n_system
                op = ADAPT.Operators.qubitexcitation(n_system, i, j)
                push!(pool, op)
                target_and_source[op] = [[i,j]]
            end
        end

        # doubles
        orbitals = collect(range(1,n_system))
        excitation_orbs = collect(combinations(orbitals,4))
        for _orbs in excitation_orbs
            target_pair = [_orbs[1],_orbs[2]]; source_pair = [_orbs[3],_orbs[4]]
            new_op = ADAPT.Operators.qubitexcitation(n_system, target_pair[1], target_pair[2], source_pair[1], source_pair[2])
            push!(pool, new_op)
            target_and_source[new_op] = [target_pair,source_pair]

            target_pair = [_orbs[1],_orbs[3]]; source_pair = [_orbs[2],_orbs[4]]
            new_op = ADAPT.Operators.qubitexcitation(n_system, target_pair[1], target_pair[2], source_pair[1], source_pair[2])
            push!(pool, new_op)
            target_and_source[new_op] = [target_pair,source_pair]

            target_pair = [_orbs[2],_orbs[3]]; source_pair = [_orbs[1],_orbs[4]]
            new_op = ADAPT.Operators.qubitexcitation(n_system, target_pair[1], target_pair[2], source_pair[1], source_pair[2])
            push!(pool, new_op)
            target_and_source[new_op] = [target_pair,source_pair]
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

    function qubitexcitationpool_complemented(n_system::Int)
        """
        Returns the complemented qubit excitation pool on n_system qubits.
        """      
        pool = ScaledPauliVector{n_system}[]
        target_and_source = Dict{ScaledPauliVector{n_system}, Vector{Vector{Int64}}}()
        for i in   1:n_system
            for j in i+1:n_system
                op = qubitexcitation_complemented(n_system, i, j)
                push!(pool, op)
                target_and_source[op] = [[i,j]]
            end
        end

        orbitals = collect(range(1,n_system))
        excitation_orbs = collect(combinations(orbitals,4))
        for _orbs in excitation_orbs
            target_pair = [_orbs[1],_orbs[2]]; source_pair = [_orbs[3],_orbs[4]]
            new_op = qubitexcitation_complemented(n_system, target_pair[1], target_pair[2], source_pair[1], source_pair[2])
            push!(pool, new_op)
            target_and_source[new_op] = [target_pair,source_pair]

            target_pair = [_orbs[1],_orbs[3]]; source_pair = [_orbs[2],_orbs[4]]
            new_op = qubitexcitation_complemented(n_system, target_pair[1], target_pair[2], source_pair[1], source_pair[2])
            push!(pool, new_op)
            target_and_source[new_op] = [target_pair,source_pair]

            target_pair = [_orbs[2],_orbs[3]]; source_pair = [_orbs[1],_orbs[4]]
            new_op = qubitexcitation_complemented(n_system, target_pair[1], target_pair[2], source_pair[1], source_pair[2])
            push!(pool, new_op)
            target_and_source[new_op] = [target_pair,source_pair]
        end
        return pool, target_and_source
    end

    function qubitadaptpool(n_system::Int)
        """
        Returns the qubit ADAPT pool on n_system qubits as defined in PRX QUANTUM 2, 020310 (2021).
        """      
        pool = ScaledPauliVector{n_system}[]
        for i in   1:n_system
            for j in i+1:n_system
                excitation_op = ADAPT.Operators.qubitexcitation(n_system, i, j)
                for sp in excitation_op
                    sp_new = ScaledPauli(1.0+0.0im, sp.pauli)
                    push!(pool, [sp_new])
                end        
            end
        end

        # doubles
        orbitals = collect(range(1,n_system))
        excitation_orbs = collect(combinations(orbitals,4))
        for _orbs in excitation_orbs
            target_pair = [_orbs[1],_orbs[2]]; source_pair = [_orbs[3],_orbs[4]]
            new_op = ADAPT.Operators.qubitexcitation(n_system, target_pair[1], target_pair[2], source_pair[1], source_pair[2])
            for sp in new_op
                sp_new = ScaledPauli(1.0+0.0im, sp.pauli)
                push!(pool, [sp_new])
            end          

            target_pair = [_orbs[1],_orbs[3]]; source_pair = [_orbs[2],_orbs[4]]
            new_op = ADAPT.Operators.qubitexcitation(n_system, target_pair[1], target_pair[2], source_pair[1], source_pair[2])
            for sp in new_op
                sp_new = ScaledPauli(1.0+0.0im, sp.pauli)
                push!(pool, [sp_new])
            end          

            target_pair = [_orbs[2],_orbs[3]]; source_pair = [_orbs[1],_orbs[4]]
            new_op = ADAPT.Operators.qubitexcitation(n_system, target_pair[1], target_pair[2], source_pair[1], source_pair[2])
            for sp in new_op
                sp_new = ScaledPauli(1.0+0.0im, sp.pauli)
                push!(pool, [sp_new])
            end          
        end
        return pool
    end

    function two_local_pool(n::Int64, axes=["I","X","Y","Z"])
        pool = ScaledPauliVector{n}[]
        for pair in Iterators.product(ntuple(i->1:n, 2)...)
            i,j = pair
            if i < j
                for pair2 in Iterators.product(ntuple(i->axes, 2)...)
                    a,b = pair2
                    if a == "I" || b == "I"  # to include 1-local strings, use: if a == b == "I"
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

    function tile_operators(L1::Int, L2::Int, chosen_operators::Array, PBCs)
        """
        Constructs the tiled operators for a system of `L2` qubits, given the operators
        chosen for a smaller problem instance on `L1` qubits.

        # Parameters
        - `L1`: number of qubits for small problem instance 
        - `L2`: number of qubits for large problem instance
        - `chosen_operators`: list of operators chosen by ADAPT for small problem instance
        - `PBCs`: periodic boundary conditions

        # Returns
        - `tiled_ops`: list of tiled operators in the form of a Vector{ScaledPauli}
        """
        @assert L2 >= L1 "L2 must be greater than L1."
        # sort the chosen operators lexically
        chosen_operators = sort(unique(chosen_operators))

        l = length(chosen_operators)
        tiled_op_strs = []
        for i=1:l
            el_str = chosen_operators[i]
            str_stripped = strip(el_str,['I'])
            permutations = PBCs ? (L2-1) : (L2-length(str_stripped))
            tiled_op_str = str_stripped*"I"^(L2-length(str_stripped))
            push!(tiled_op_strs,tiled_op_str)
            for j = 1:permutations
                new_s = tiled_op_str[end]*tiled_op_str[1:end-1]
                tiled_op_str = new_s
                push!(tiled_op_strs,tiled_op_str)     
            end
        end
        tiled_op_strs = sort(unique(tiled_op_strs))

    #     println(tiled_op_strs)
        tiled_ops = ScaledPauliVector{L2}[]
        for str in tiled_op_strs
            push!(tiled_ops,[ScaledPauli(Pauli(str))])
        end

        return tiled_ops
    end

    function tile_op_strings(L1::Int, L2::Int, chosen_operators::Array, PBCs)
        """
        Constructs the tiled operators for a system of `L2` qubits, given the operators
        chosen for a smaller problem instance on `L1` qubits.

        # Parameters
        - `L1`: number of qubits for small problem instance 
        - `L2`: number of qubits for large problem instance
        - `chosen_operators`: list of operators chosen by ADAPT for small problem instance
        - `PBCs`: periodic boundary conditions

        # Returns
        - `tiled_ops`: list of tiled operators in the form of a Vector{String}
        """
        @assert L2 >= L1 "L2 must be greater than L1."
        # sort the chosen operators lexically
        chosen_operators = sort(unique(chosen_operators))

        l = length(chosen_operators)
        tiled_op_strs = []
        for i=1:l
            el_str = chosen_operators[i]
            str_stripped = strip(el_str,['I'])
            permutations = PBCs ? (L2-1) : (L2-length(str_stripped))
            tiled_op_str = str_stripped*"I"^(L2-length(str_stripped))
            push!(tiled_op_strs,tiled_op_str)
            for j = 1:permutations
                new_s = tiled_op_str[end]*tiled_op_str[1:end-1]
                tiled_op_str = new_s
                push!(tiled_op_strs,tiled_op_str)     
            end
        end
        tiled_op_strs = sort(unique(tiled_op_strs))
        return tiled_op_strs
    end

    function tile_ops(L1::Int, L2::Int, chosen_operators::Vector{Vector{ScaledPauli{N}}}, PBCs) where {N}
        """
        Constructs the tiled operators for a system of `L2` qubits, given the operators
        chosen for a smaller problem instance on `L1` qubits.

        # Parameters
        - `L1`: number of qubits for small problem instance 
        - `L2`: number of qubits for large problem instance
        - `chosen_operators`: list of operators chosen by ADAPT for small problem instance
        - `PBCs`: periodic boundary conditions

        # Returns
        - `tiled_ops`: list of tiled operators in the form of a Vector{ScaledPauli}
        """
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
        # println(tiled_ops)
        return tiled_ops
    end
end