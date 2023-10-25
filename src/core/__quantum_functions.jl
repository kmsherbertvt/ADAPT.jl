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
    evolve_state!(
        ansatz::AbstractAnsatz,
        state::QuantumState,
    )

Apply an ansatz to the given quantum state, mutating the state.

"""
function evolve_state!(
    ansatz::AbstractAnsatz,
    state::QuantumState,
)::QuantumState
    generators = get_generators(ansatz)
    parameters = get_parameters(ansatz)
    for i in eachindex(generators)
        evolve_state!(generators[i], parameters[i], state)
    end
    return state # It's conventional to return the mutated object, but it's not mandatory.
end

"""
    evolve_state!(
        G::Generator,
        θ::Parameter,
        ψ::QuantumState,
    )

Rotate a quantum state `ψ` by an amount `x` about the axis defined by `G`, mutating `ψ`.

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
    calculate_gradient(
        ansatz::AbstractAnsatz,
        H::Observable,
        ψ0::QuantumState,
    )

Evaluate the gradient of a cost-function with respect to a particular ADAPT state.

# Parameters
- `ansatz`: the ADAPT state
- `H`: the object defining the cost-function
- `ψ0`: an initial quantum state which the `ansatz` operates on

# Returns
- `gradient`: a vector of type `typeof_energy(observable)`.
        `gradient[i]` is `∂⟨H⟩/∂x[i]`, where `x` is the vector `get_parameters(ansatz)`

"""
function calculate_gradient(
    ansatz::AbstractAnsatz,
    observable::Observable,
    reference::QuantumState,
)
    L = length(get_parameters(ansatz))
    gradient = Vector{energy_type(observable)}(undef, L)
    calculate_gradient!(gradient, ansatz, observable, reference)
    return gradient
end


"""
    calculate_gradient!(
        gradient::EnergyList,
        ansatz::AbstractAnsatz,
        observable::Observable,
        reference::QuantumState,
    )

Update `gradient` to that of a cost-function with respect to a particular ADAPT state.

# Implementation

This is simply an "in-place" version of `calculate_gradient`;
    the idea is that you only have to allocate the gradient vector once,
    at the start of an optimization.

Most high-performance optimization packages prefer this sort of interface,
    so it is the one you should define, and it is usually the one you should use.

The "not-in-place" version of the method is really
    only meant for convenience when playing in the REPL.

"""
calculate_gradient!(
    gradient::EnergyList,
    ansatz::AbstractAnsatz,
    observable::Observable,
    reference::QuantumState,
) = NotImplementedError()