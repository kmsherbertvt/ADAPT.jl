#= DISCLAIMER: This file is meant to be included in ADAPT.jl.
    Names not defined in this file or in base Julia can be found there.
=#

import PauliOperators

"""
    Generator

The Union of any type that could be used as a pool operator.

# Implemented Types

Any type at all can be used as a generator
  if there is a compatible implementation
  of the methods listed in the `Implementation` section.

The following types have implementations fleshed out in this library already:
- `PauliOperators.AbstractPauli`: A single Pauli word
- `PauliOperators.PauliSum`: A Hermitian operator decomposed into the Pauli basis
- `PauliOperators.ScaledPauliVector`: Same but with a different internal data structure

  For each of the above,
      the generator `G` generates the unitary `exp(-iθG)`.
  Hermiticity in `G` is not enforced, so be careful when constructing your pool operators.

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
    PauliOperators.AbstractPauli,
    PauliOperators.PauliSum,
    PauliOperators.ScaledPauliVector,
    Any,
}

"""
    GeneratorList

Semantic alias for a vector of generators.

"""
GeneratorList = AbstractVector{<:Generator}



"""
    Observable

The Union of any type that could define a cost function.

The type is so named because the typical cost-function is the expecation value
    of a Hermitian operator, aka a quantum observable.

# Implemented Types

Any type at all can be used as a generator
  if there is a compatible implementation
  of the methods listed in the `Implementation` section.

The following types have implementations fleshed out in this library already:
- `PauliOperators.AbstractPauli`: A single Pauli word
- `PauliOperators.PauliSum`: A Hermitian operator decomposed into the Pauli basis
- `PauliOperators.ScaledPauliVector`: Same but with a different internal data structure

  For each of the above,
      the evaluation of `H` with respect to a quantum state `|Ψ⟩` is `⟨Ψ|H|Ψ⟩`.

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
    PauliOperators.AbstractPauli,
    PauliOperators.PauliSum,
    PauliOperators.ScaledPauliVector,
    Any,
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
    QuantumState

The Union of any type that could define a quantum state.

# Implemented Types

Any type at all can be used as a generator
  if there is a compatible implementation
  of the methods listed in the `Implementation` section.

The following types have implementations fleshed out in this library already:
- `Vector{<:Complex}`: A dense statevector in the computational basis
- `PauliOperators.SparseKetBasis`:

  A dict mapping individual kets (`PauliOperators.KetBitString`) to their coefficients.

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
    PauliOperators.SparseKetBasis,
    Vector{<:Complex},
    Any,
}
