#= DISCLAIMER: This file is meant to be included in ADAPT.jl.
    Names not defined in this file or in base Julia can be found there.
=#

"""
    AbstractAnsatz{F,G}

An image of an ADAPT protocol in a frozen state.

The type is so named because the most basic possible such image
    consists of just the generators and parameters of the ansatz thus far selected,
    but richer variants of ADAPT will have richer states.

For example, a version of ADAPT which carries information
    on the inverse Hessian across ADAPT iterations
    would need to operate on an ansatz type which includes the inverse Hessian.

Nevertheless, every sub-type of AbstractAnsatz implements the AbstractVector interface,
    where elements are pairs `(generator => parameter)`.

So, for example, an ansatz maintaining an inverse Hessian would need
    to override `push!` `insert!`, etc. to ensure the dimension of the Hessian matches.

# Type Parameters
- `F`: the number type of the parameters (usually `Float64`)
- `G`: the subtype of `Generator`

# Implementation

Sub-types must implement the following methods:
- `__get__generators(::AbstractAnsatz{F,G})::Vector{G}`
- `__get__parameters(::AbstractAnsatz{F,G})::Vector{F}`
- `__get__optimized(::AbstractAnsatz{F,G})::Ref{Bool}`
- `__get__converged(::AbstractAnsatz{F,G})::Ref{Bool}`

Each of these is expected to simply retrieve an attribute of the struct.
You can call them whatever you'd like, but functionally, here's what they mean:
- `generators::Vector{G}`: the sequence of generators
- `parameters::Vector{F}`: the corresponding sequence of parameters

    Note that these vectors will be mutated and resized as needed.

- `optimized::Ref{Bool}`: a flag indicating that the current parameters are optimal
- `converged::Ref{Bool}`: a flag indicating that the current generators are optimal

    Note that these must be of type `Ref`, so their values can be toggled as needed.

In addition, there must be a compatible implementation for each of:

- `partial(
    index::Int,
    ansatz::AbstractAnsatz,
    observable::Observable,
    reference::QuantumState,
  )`

- `calculate_score(
    ::AbstractAnsatz,
    ::AdaptProtocol,
    ::Generator,
    ::Observable,
    ::QuantumState,
  )::Score`

- `adapt!(
    ::AbstractAnsatz,
    ::Trace,
    ::AdaptProtocol,
    ::GeneratorList,
    ::Observable,
    ::QuantumState,
    ::CallbackList,
  )`

- `optimize!(
    ::AbstractAnsatz,
    ::Trace,
    ::OptimizationProtocol,
    ::Observable,
    ::QuantumState,
    ::CallbackList,
  )`

That said, most basic implementations of these methods are defined for abstract ansatze,
    so you oughtn't need to worry about them.

Please see individual method documentation for details.

"""
abstract type AbstractAnsatz{F<:Parameter,G<:Generator} <: AbstractVector{Pair{G,F}} end

__get__generators(::AbstractAnsatz) = NotImplementedError()
__get__parameters(::AbstractAnsatz) = NotImplementedError()
__get__optimized(::AbstractAnsatz) = NotImplementedError()
__get__converged(::AbstractAnsatz) = NotImplementedError()


##########################################################################################
#= AbstractVector interface. =#

Base.size(ansatz::AbstractAnsatz) = size(__get__generators(ansatz))
Base.IndexStyle(::Type{<:AbstractAnsatz}) = IndexLinear()

function Base.getindex(ansatz::AbstractAnsatz, i::Int)
    return __get__generators(ansatz)[i] => __get__parameters(ansatz)[i]
end

function Base.setindex!(ansatz::AbstractAnsatz{F,G}, pair::Pair{G,F}, i::Int) where {F,G}
    __get__generators(ansatz)[i] = pair.first
    __get__parameters(ansatz)[i] = pair.second
end

function Base.resize!(ansatz::AbstractAnsatz, nl::Int)
    resize!(__get__generators(ansatz), nl)
    resize!(__get__parameters(ansatz), nl)
end

##########################################################################################
#= Ansatz-specific interface. =#

"""
    typeof_parameter(::AbstractAnsatz)

The number type of the variational parameters in this ansatz.

I think this will always be a sub-type of `AbstractFloat`,
    and almost always `Float64`.

"""
typeof_parameter(::AbstractAnsatz{F,G}) where {F,G} = F

"""
    is_optimized(::AbstractAnsatz)

Check whether the ansatz parameters are flagged as optimal.

Note that this is a state variable in its own right;
    its value is independent of the actual parameters themselves,
    but depends on all the protocols and callbacks
    which brought the ansatz to its current state.

"""
is_optimized(ansatz::AbstractAnsatz) = __get__optimized(ansatz)[]


"""
    is_converged(::AbstractAnsatz)

Check whether the sequence of generators in this ansatz are flagged as optimal.

Note that this is a state variable in its own right;
    its value is independent of the actual generators themselves,
    but depends on all the protocols and callbacks
    which brought the ansatz to its current state.

"""
is_converged(ansatz::AbstractAnsatz) = __get__converged(ansatz)[]

"""
    set_optimized!(::AbstractAnsatz, ::Bool)

Flag the ansatz parameters as optimal.

"""
set_optimized!(ansatz::AbstractAnsatz, flag::Bool) = __get__optimized(ansatz)[] = flag


"""
    set_converged!(::AbstractAnsatz, ::Bool)

Flag the sequence of generators in this ansatz as optimal.

"""
set_converged!(ansatz::AbstractAnsatz, flag::Bool) = __get__converged(ansatz)[] = flag

"""
    angles(::AbstractAnsatz)

Fetch all parameters in the ansatz as a vector.

"""
angles(ansatz::AbstractAnsatz) = __get__parameters(ansatz)

"""
    bind!(::AbstractAnsatz, ::ParameterList)

Replace all parameters in the ansatz.

"""
function bind!(ansatz::AbstractAnsatz{F,G}, x::AbstractVector{F}) where {F,G}
    __get__parameters(ansatz) .= x
end