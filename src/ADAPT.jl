module ADAPT
    #= TODO: Package and manifest are over-burdened.
            Ask Nick to re-organize the dependencies in PauliOperators.jl.
    =#

    # A stock error for semantic implementation of abstract methods.
    NotImplementedError() = error("Method not yet implemented for these arguments.")

    # Alias `Number` a few times to get semantic parameter types.
    include("core/__number_aliases.jl")
    export Parameter, ParameterList
    export Energy, EnergyList
    export Score, ScoreList

    # Define types to represent quantum objects (operators and states).
    include("core/__quantum_objects.jl")
    export Generator, GeneratorList
    export Observable, typeof_energy
    export QuantumState

    # Define the type for the ADAPT state (which I call an "ansatz").
    include("core/__ansatz.jl")
    export AbstractAnsatz
    export typeof_parameter
    export is_optimized, set_optimized!
    export is_converged, set_converged!
    export angles, bind!

    # Define functions for quantum operations (evolution and expectation).
    include("core/__quantum_functions.jl")
    export evolve_state, evolve_state!
    export evaluate, partial, gradient, gradient!
    export make_costfunction, make_gradfunction, make_gradfunction!

    # Define functions for dense matrix representations.
    include("core/__matrix_functions.jl")
    export evolve_unitary, evolve_unitary!

    # Define the interface for implementing distinct ADAPT and VQE protocols.
    include("core/__protocols.jl")
    export Trace, Data
    export AbstractCallback, CallbackList
    export AdaptProtocol, adapt!
    export typeof_score, calculate_score, calculate_scores
    export OptimizationProtocol, optimize!
    export run!

    # Define standard functions to validate user-defined types work with ADAPT.
    include("core/__validation.jl")
    export validate
    export validate_runtime, validate_consistency
    export validate_evolution, validate_evaluation, validate_gradient, validate_scores

    module Basics
        # A minimalist's ansatz. Sufficient for most applications.
        include("base/Ansatz.jl")
        export Ansatz

        # The original ADAPT: scoring by expectation value of commutators
        include("base/VanillaADAPT.jl")
        export VANILLA

        # Interface into any optimization method implemented by `Optim.jl`
        include("base/OptimOptimizer.jl")
        export OptimOptimizer

        # Stuff that should be in PauliOperators.jl
        module MyPauliOperators; include("base/__paulioperators.jl"); end
        # Implementing generators and observables with our group's `PauliOperators.jl`
        include("base/pauli_plugin.jl")

        # A suite of callbacks for tracing, printing, and stopping
        include("base/Callbacks.jl")
        export Callbacks

        # A suite of common operators, especially useful for constructing operator pools.
        include("base/Operators.jl")
        export Operators

        # A suite of common lattice Hamiltonians.
        include("hamiltonians/latticemodels.jl")
        export LatticeModelHamiltonians

        # A suite of common operator pools.
        include("pools/pools.jl")
        export OperatorPools
    end
    using .Basics
    export Ansatz
    export VANILLA
    export OptimOptimizer
    export Callbacks
    export Operators
    export LatticeModelHamiltonians
    export OperatorPools

    module OptimizationFreeADAPT
        include("optimizationfree/OptimizationFree.jl")
        export OPTIMIZATION_FREE
    end


    module OverlapADAPT
        include("overlap/OverlapADAPT.jl")
        export Infidelity
    end
end

#=

TODO:
- Some pools. Grab from Diksha's code.
- Lie rank calculation, in pauli-dedicated code.

=#