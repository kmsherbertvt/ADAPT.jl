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
    export AbstractGenerator, Generator, GeneratorList
    export AbstractObservable, Observable, typeof_energy
    export AbstractQuantumState, QuantumState

    # Define the type for the ADAPT state (which I call an "ansatz").
    include("core/__ansatz.jl")
    export AbstractAnsatz
    export typeof_parameter
    export is_optimized!, set_optimized!
    export is_converged, set_converged!
    export update_parameters!

    # Define functions for quantum operations (evolution and expectation).
    include("core/__quantum_functions.jl")
    export evolve_state, evolve_state!
    export evaluate, partial

    # Define the interface for implementing distinct ADAPT and VQE protocols.
    include("core/__protocols.jl")
    export Trace, Data
    export AbstractCallback, CallbackList
    export AdaptProtocol, adapt!
    export typeof_score, calculate_score, calculate_scores
    export OptimizationProtocol, optimize!
    export run!

    module Basics
        # Implementing generators and observables with our group's `PauliOperators.jl`
        include("base/pauli_plugin.jl")

        # A minimalist's ansatz. Sufficient for most applications.
        include("base/Ansatz.jl")
        export Ansatz

        # The original ADAPT: scoring by expectation value of commutators
        include("base/VanillaADAPT.jl")
        export VANILLA

        # Interface into any optimization method implemented by `Optim.jl`
        include("base/OptimOptimizer.jl")
        export OptimOptimizer

        # A minimalist's optimization protocol.
        include("base/OptimizationFree.jl")
        export OPTIMIZATION_FREE

        # A suite of callbacks for tracing, printing, and stopping
        include("base/Callbacks.jl")
        export Callbacks
    end
    using .Basics
    export Ansatz
    export VANILLA
    export OptimOptimizer, OPTIMIZATION_FREE
    export Callbacks
end