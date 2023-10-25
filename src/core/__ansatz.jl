#= DISCLAIMER: This file is meant to be included in ADAPT.jl.
    Names not defined in this file or in base Julia can be found there.
=#

"""
    AbstractAnsatz

An image of an ADAPT protocol in a frozen state.

The type is so named because the most basic possible such image
    consists of just the generators and parameters of the ansatz thus far selected,
    but richer variants of ADAPT will have richer states.

For example, a version of ADAPT which carries information
    on the inverse Hessian across ADAPT iterations
    would need to operate on an ansatz type which includes the inverse Hessian.

# Implementation

Sub-types must implement the following methods:
- `typeof_parameter(::AbstractAnsatz)::Type{<:Parameter}`
- `get_generators(::AbstractAnsatz)::GeneratorList`
- `get_parameters(::AbstractAnsatz)::ParameterList`
- `is_optimized(::AbstractAnsatz)::Bool`
- `is_converged(::AbstractAnsatz)::Bool`
- `set_optimized!(::AbstractAnsatz, ::Bool)`
- `set_converged!(::AbstractAnsatz, ::Bool)`
- `add_generator!(::AbstractAnsatz, ::Generator, ::Parameter)`
- `update_parameters!(::AbstractAnsatz, ::ParameterList)`

In addition, there must be a compatible implementation for each of:
- `calculate_gradient!(
    ::EnergyList,
    ::AbstractAnsatz,
    ::Observable,
    ::QuantumState,
  )`

- `calculate_score(
    ::AbstractAnsatz,
    ::AdaptProtocol,
    ::Generator,
    ::Observable,
    ::QuantumState,
  )::Score`

- `optimize!(
    ::AbstractAnsatz,
    ::Trace,
    ::OptimizationProtocol,
    ::Observable,
    ::QuantumState,
    ::CallbackList,
  )`

Please see individual method documentation for details.

"""
abstract type AbstractAnsatz end

"""
    typeof_parameter(::AbstractAnsatz)

The number type of the variational parameters in this ansatz.

I think this will always be a sub-type of `AbstractFloat`,
    and almost always `Float64`.

"""
typeof_parameter(::AbstractAnsatz) = NotImplementedError()

"""
    get_parameters(::AbstractAnsatz)

The vector of parameters currently used in this ansatz.

# Implementation

Return a vector whose elements sub-type `Parameter`,
    presumably just fetching an attribute of your object.

To allow protocols and callbacks the greatest flexibility,
    modifications to the vector returned here SHOULD modify the ansatz itself.

"""
get_parameters(::AbstractAnsatz) = NotImplementedError()

"""
    get_generators(::AbstractAnsatz)

The vector of parameters currently used in this ansatz.

# Implementation

Return a vector whose elements sub-type `Generator`,
    presumably just fetching an attribute of your object.

To allow protocols and callbacks the greatest flexibility,
    modifications to the vector returned here SHOULD modify the ansatz itself.

"""
get_generators(::AbstractAnsatz) = NotImplementedError()

"""
    is_optimized(::AbstractAnsatz)

Check whether the ansatz parameters are flagged as optimal.

Note that this is a state variable in its own right;
    its value is independent of the actual parameters themselves,
    but depends on all the protocols and callbacks
    which brought the ansatz to its current state.

# Implementation

Return a Bool, presumably just fetching an attribute of your type.

"""
is_optimized(::AbstractAnsatz) = NotImplementedError()

"""
    is_converged(::AbstractAnsatz)

Check whether the sequence of generators in this ansatz are flagged as converged.

Note that this is a state variable in its own right;
    its value is independent of the actual generators themselves,
    but depends on all the protocols and callbacks
    which brought the ansatz to its current state.

# Implementation

Return a Bool, presumably just fetching an attribute of your type.

"""
is_converged(::AbstractAnsatz) = NotImplementedError()

"""
    set_optimized!(::AbstractAnsatz, ::Bool)

Flag the ansatz parameters as optimal.

# Implementation

Presumably, just set an attribute of your type.

"""
set_optimized!(::AbstractAnsatz, ::Bool) = NotImplementedError()

"""
    set_converged!(::AbstractAnsatz, ::Bool)

Flag the sequence of generators in this ansatz as converged.

# Implementation

Presumably, just set an attribute of your type.

"""
set_converged!(::AbstractAnsatz, ::Bool) = NotImplementedError()

"""
    add_generator!(::AbstractAnsatz, ::Generator, ::Parameter)

Add a new generator (with corresponding parameter) to the ansatz.

# Implementation

Presumably, just append each argument to the appropriate vector attribute of your type.

"""
add_generator!(::AbstractAnsatz, ::Generator, ::Parameter) = NotImplementedError()

"""
    update_parameters!(::AbstractAnsatz, ::ParameterList)

Replace all parameters in the ansatz.

# Implementation

Presumably, just copy the given parameter vector into a vector attribute of your type.

"""
update_parameters!(::AbstractAnsatz, ::ParameterList) = NotImplementedError()

