#= DISCLAIMER: This file is meant to be included in ADAPT.jl.
    Names not defined in this file or in base Julia can be found there.
=#

"""
    evolve_state(
        ansatz::AbstractAnsatz,
        reference::QuantumState,
    )

Calculate the quantum state resulting from applying an ansatz to a given reference state.

"""
function evolve_state(
    ansatz::AbstractAnsatz,
    reference::QuantumState,
)
    state = deepcopy(reference)
    evolve_state!(ansatz, state)
    return state
end

"""
    evolve_state(
        G::Generator,
        θ::Parameter,
        ψ::QuantumState,
    )

Calculate the quantum state rotating `ψ` by an amount `x` about the axis defined by `G`.

"""
function evolve_state(
    G::Generator,
    θ::Parameter,
    Ψ::QuantumState,
)
    state = deepcopy(Ψ)
    evolve_state!(G, θ, state)
    return state
end

"""
    evolve_state!(
        ansatz::AbstractAnsatz,
        state::QuantumState,
    )

Apply an ansatz to the given quantum state, mutating and returning the state.

By default, generators with a lower index are applied to the state earlier.
This means that the equation for |Ψ⟩ would list out generators in reverse order.
Specific implementations of ansatze may override this behavior.

"""
function evolve_state!(
    ansatz::AbstractAnsatz,
    state::QuantumState,
)::QuantumState
    for (generator, parameter) in ansatz
        evolve_state!(generator, parameter, state)
    end
    return state
end

"""
    evolve_state!(
        G::Generator,
        θ::Parameter,
        ψ::QuantumState,
    )

Rotate a quantum state `ψ` by an amount `x` about the axis defined by `G`,
    mutating and returning `ψ`.

# Implementation

Typically, the "rotation" is the unitary operator `exp(-iθG)`,
    but different `Generator` types could have different effects.

"""
evolve_state!(
    ::Generator,
    ::Parameter,
    ::QuantumState,
) = NotImplementedError()




"""
    evaluate(
        ansatz::AbstractAnsatz,
        H::Observable,
        ψ0::QuantumState,
    )

Evaluate a cost-function with respect to a particular ADAPT state.

# Parameters
- `ansatz`: the ADAPT state
- `H`: the object defining the cost-function
- `ψ0`: an initial quantum state which the `ansatz` operates on

# Returns
- `energy`: a scalar number, whose type is `typeof_energy(observable)`

"""
function evaluate(
    ansatz::AbstractAnsatz,
    observable::Observable,
    reference::QuantumState,
)
    state = evolve_state(ansatz, reference)
    return evaluate(observable, state)
end

"""
    evaluate(
        H::Observable,
        Ψ::QuantumState,
    )

Evaluate a cost-function with respect to a particular quantum state.

# Parameters
- `H`: the object defining the cost-function
- `Ψ`: the quantum state

# Returns
- `energy`: a scalar number, whose type is `typeof_energy(observable)`

# Implementation

Typically, the "cost-function" is the expectation value `⟨Ψ|H|Ψ⟩`,
    but different `Observable` types could have different definitions.

"""
evaluate(
    ::Observable,
    ::QuantumState,
) = NotImplementedError()


"""
    partial(
        index::Int,
        ansatz::AbstractAnsatz,
        observable::Observable,
        reference::QuantumState,
    )

The partial derivative of a cost-function with respect to the i-th parameter in an ansatz.

# Parameters
- `index`: the index of the parameter to calculate within `ansatz`
- `ansatz`: the ADAPT state
- `H`: the object defining the cost-function
- `ψ0`: an initial quantum state which the `ansatz` operates on

# Returns
- a number of type `typeof_energy(observable)`.

# Implementation

Typically, generators apply a unitary rotation,
    so the partial consists of a partial evolution up to the indexed generator,
    then a "kick" from the generator itself, then a final evolution,
    and a braket with the observable.
But, different ansatze may have a different procedure.

"""
partial(
    index::Int,
    ansatz::AbstractAnsatz,
    observable::Observable,
    reference::QuantumState,
) = NotImplementedError()

"""
    gradient(
        ansatz::AbstractAnsatz,
        observable::Observable,
        reference::QuantumState,
    )

Construct a vector of partial derivatives with respect to each parameter in the ansatz.

# Parameters
- `ansatz`: the ADAPT state
- `H`: the object defining the cost-function
- `ψ0`: an initial quantum state which the `ansatz` operates on

# Returns
- a vector of type `typeof_energy(observable)`.

"""
function gradient(
    ansatz::AbstractAnsatz,
    observable::Observable,
    reference::QuantumState,
)
    F = typeof_energy(observable)
    result = Vector{F}(undef, length(ansatz))
    result = gradient!(result, ansatz, observable, reference)
    return result
end

"""
    gradient!(
        result::AbstractVector,
        ansatz::AbstractAnsatz,
        observable::Observable,
        reference::QuantumState,
    )

Fill a vector of partial derivatives with respect to each parameter in the ansatz.

# Parameters
- `result`: vector which will contain the gradient after calling this function
- `ansatz`: the ADAPT state
- `H`: the object defining the cost-function
- `ψ0`: an initial quantum state which the `ansatz` operates on

# Returns
- `result`

"""
function gradient!(
    result::AbstractVector,
    ansatz::AbstractAnsatz,
    observable::Observable,
    reference::QuantumState,
)
    for i in eachindex(ansatz)
        result[i] = ADAPT.partial(i, ansatz, observable, reference)
    end
    return result
end