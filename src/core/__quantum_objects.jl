#= DISCLAIMER: This file is meant to be included in ADAPT.jl.
    Names not defined in this file or in base Julia can be found there.
=#

import PauliOperators: AbstractPauli
import PauliOperators: SparseKetBasis

"""
    AbstractGenerator

An abstract type included in the `Generator` Union,
    so that users may implement their own without having to modify the package.

"""
abstract type AbstractGenerator end

"""
    Generator

The Union of any type that could be used as a pool operator.

# Implementation

This list should be extended as new sub-types are implemented.

There must be a compatible implementation for each of:

- `evolve_state!(
    ::Generator,
    ::Parameter,
    ::QuantumState,
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

"""
Generator = Union{
    AbstractGenerator,
    AbstractPauli,
}

"""
    GeneratorList

Semantic alias for a vector of generators.

"""
GeneratorList = AbstractVector{<:Generator}


"""
    AbstractObservable

An abstract type included in the `Observable` Union,
    so that users may implement their own without having to modify the package.

"""
abstract type AbstractObservable end

"""
    Observable

The Union of any type that could define a cost function.

The type is so named because the typical cost-function is the expecation value
    of a Hermitian operator, aka a quantum observable.

# Implementation

This list should be extended as new sub-types are implemented.

Sub-types must implement the following method:
- `typeof_energy(::Observable)::Type{<:Energy}`

In addition, there must be a compatible implementation for each of:

- `evaluate(
    ::Observable,
    ::QuantumState,
  )::Energy`

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
  )::Bool`

"""
Observable = Union{
    AbstractObservable,
    AbstractPauli,
}

"""
    typeof_energy(::Observable)

The number type of a cost-function.

The method is so named because the typical cost-function is the expecation value
    of a Hamiltonian, aka an energy.

# Implementation

Usually a sub-type of `AbstractFloat`, and probably just about always `Float64`.

"""
typeof_energy(::Observable) = NotImplementedError()


"""
    AbstractQuantumState

An abstract type included in the `QuantumState` Union,
    so that users may implement their own without having to modify the package.

"""
abstract type AbstractQuantumState end

"""
    QuantumState

The Union of any type that could define a quantum state.

# Implementation

This list should be extended as new sub-types are implemented.

There must be a compatible implementation for each of:

- `evolve_state!(
    ::Generator,
    ::Parameter,
    ::QuantumState,
  )`

- `evaluate(
    ::Observable,
    ::QuantumState,
  )::Energy`

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
  )::Bool`

"""
QuantumState = Union{
    AbstractQuantumState,
    SparseKetBasis,
    Vector{<:Complex},
}
