module ADAPT
    NotImplementedError() = error("Method not yet implemented for these arguments.")

    include("core/__number_aliases.jl")
    export Parameter, ParameterList
    export Energy, EnergyList
    export Score, ScoreList

    include("core/__quantum_objects.jl")
    export AbstractGenerator, Generator, GeneratorList
    export AbstractObservable, Observable, typeof_energy
    export AbstractQuantumState, QuantumState

    include("core/__ansatz.jl")
    export AbstractAnsatz
    export typeof_parameter
    export get_generators, get_parameters
    export is_optimized!, set_optimized!
    export is_converged, set_converged!
    export add_generator!, update_parameters!

    include("core/__quantum_functions.jl")
    export evolve_state, evolve_state!
    export evaluate
    export calculate_gradient, calculate_gradient!

    include("core/__callbacks.jl")
    export Trace, Data
    export AbstractCallback, CallbackList

    include("core/__protocols.jl")
    export AdaptProtocol, adapt!
    export typeof_score, calculate_score, calculate_scores
    export OptimizationProtocol, optimize!
    export run!



    module Base
        include("base/Ansatz.jl")
        export Ansatz

        include("base/OptimizationProtocols.jl")
        export OPTIMIZATION_FREE

    end
    using Base
    export Ansatz
    export OPTIMIZATION_FREE


end