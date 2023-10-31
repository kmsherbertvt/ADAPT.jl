import ..ADAPT

"""
    Ansatz{F<:Parameter,G<:Generator}(
        parameters::Vector{F},
        generators::Vector{G},
        optimized::Bool,
        converged::Bool,
    )

A minimal ADAPT state.

# Type Parameters
- `F`: the number type for the parameters (usually `Float64` is appropriate.)
- `G`: the generator type. Any type will do, but it's best to be specific.

# Parameter
- `parameters`: list of current parameters
- `generators`: list of current generators
- `optimized`: whether the current parameters are flagged as optimal
- `converged`: whether the current generators are flagged as converged

"""
struct Ansatz{F,G} <: ADAPT.AbstractAnsatz{F,G}
    generators::Vector{G}
    parameters::Vector{F}
    optimized::Ref{Bool}
    converged::Ref{Bool}
end

ADAPT.__get__generators(ansatz::Ansatz) = ansatz.generators
ADAPT.__get__parameters(ansatz::Ansatz) = ansatz.parameters
ADAPT.__get__optimized(ansatz::Ansatz) = ansatz.optimized
ADAPT.__get__converged(ansatz::Ansatz) = ansatz.converged

"""
    Ansatz(F, G)

Convenience constructor for initializing an empty ansatz.

# Parameters
- the parameter type OR an instance of that type OR a vector whose elements are that type
- the generator type OR an instance of that type OR a vector whose elements are that type

The easiest way to use this constructor is probably to prepare your generator pool first,
    then call `Ansatz(Float64, pool)`.
But please note, the ansatz is always initialized as empty,
    even though you've passed a list of generators in the constructor!

"""
Ansatz(F, G) = Ansatz(eltype(G)[], eltype(F)[], Ref(true), Ref(false))