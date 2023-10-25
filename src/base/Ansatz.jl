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
struct Ansatz{F,G} <: ADAPT.AbstractAnsatz where {
    F <: ADAPT.Parameter,
    G <: ADAPT.Generator,
}
    parameters::Vector{F}
    generators::Vector{G}
    optimized::Bool
    converged::Bool
end

"""
    Ansatz(::F, ::G)

Convenience constructor for initializing an empty ansatz

# Parameters
- the parameter type OR an instance of that type OR a vector whose elements are that type
- the generator type OR an instance of that type OR a vector whose elements are that type

The easiest way to use this constructor is probably to prepare your generator pool first,
    then call `Ansatz(Float64, pool)`.
But please note, the ansatz is always initialized as empty,
    even though you've passed a list of generators in the constructor!

"""
Ansatz(::F, ::G) where {F,G} = Ansatz{Float64}(eltype(F)[], eltype(G)[], true, false)

ADAPT.typeof_parameter(::Ansatz{F,G}) where {F,G} = F
ADAPT.get_parameters(ansatz::Ansatz) = ansatz.parameters
ADAPT.get_generators(ansatz::Ansatz) = ansatz.generators
ADAPT.is_optimized(ansatz::Ansatz) = ansatz.optimized
ADAPT.is_converged(ansatz::Ansatz) = ansatz.converged
ADAPT.set_optimized!(ansatz::Ansatz, flag::Bool) = (ansatz.optimized = flag)
ADAPT.set_converged!(ansatz::Ansatz, flag::Bool) = (ansatz.converged = flag)

function ADAPT.add_generator!(ansatz::Ansatz{F,G}, generator::G, parameter::F) where {F,G}
    push!(ansatz.parameters, parameter)
    push!(ansatz.generators, generator)
end

function ADAPT.update_parameters!(ansatz::Ansatz{F,G}, parameters::AbstractVector{F})
    ansatz.parameters .= parameters
end