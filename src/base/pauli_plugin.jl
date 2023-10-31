import ..ADAPT

import .MyPauliOperators
import .MyPauliOperators: SparseKetBasis
import .MyPauliOperators: AbstractPauli, FixedPhasePauli, ScaledPauliVector, PauliSum
PauliOperators = MyPauliOperators
# TODO: Replace `MyPauliOperators` with `PauliOperators`

AnyPauli = Union{AbstractPauli, PauliSum, ScaledPauliVector}

ADAPT.typeof_energy(::AnyPauli) = Float64
#= NOTE:
    This type assertion assumes two things:
    1. Scaled AbstractPaulis are hard-coded to used `ComplexF64`.
    2. The observable is Hermitian (so its energy is guaranteed to be real).

    Strictly speaking, using the plain `Pauli` type undermines both assumptions,
        since half of them are anti-Hermitian and they don't have *any* float type.
    So...try to avoid using `Pauli` as an `Observable`?
    It should work fine if the phase is real, but...it'd be less than robust...
=#

function ADAPT.evolve_state!(
    G::FixedPhasePauli,
    θ::ADAPT.Parameter,
    Ψ::ADAPT.QuantumState,
)
    PauliOperators.rotate!(Ψ, G, θ)
    Ψ isa SparseKetBasis && PauliOperators.clip!(Ψ)
end

function ADAPT.evolve_state!(
    G::ScaledPauliVector,
    θ::ADAPT.Parameter,
    Ψ::ADAPT.QuantumState,
)
    for scaled in G
        # CALCULATE COEFFICIENT, ACCOUNTING FOR IMPLICIT `i` FOR EACH Y in `G.pauli`
        coeff = PauliOperators.get_phase(scaled.pauli) * scaled.coeff
        # EVOLVE BY FIXED-PHASE PAULI
        ADAPT.evolve_state!(scaled.pauli, coeff * θ, Ψ)
    end
end

function ADAPT.evolve_state!(
    G::PauliSum,
    θ::ADAPT.Parameter,
    Ψ::ADAPT.QuantumState,
)
    #= TODO: Implement this, using KrylovKit.exponentiate!

    https://jutho.github.io/KrylovKit.jl/stable/man/matfun/#KrylovKit.exponentiate

    This requires `x` to be "any Julia type with vector like behavior".
    Not too sure *which* behavior is required.
    So, there's a chance `PauliSum` will only be compatible with `AbstractVector`,
        but that's okay.

    =#
    return NotImplementedError()
end



function ADAPT.evaluate(
    H::AnyPauli,
    Ψ::ADAPT.QuantumState,
)
    return real(PauliOperators.expectation_value(H, Ψ))
end


"""
    partial(
        index::Int,
        ansatz::AbstractAnsatz,
        observable::Observable,
        reference::QuantumState,
    )

The partial derivative of a cost-function with respect to the i-th parameter in an ansatz.

The ansatz is assumed to apply a unitary rotation `exp(-iθG)`,
    where `G` is the (Hermitian) generator,
    and generators with a lower index are applied to the state earlier.
Ansatz sub-types may change both behaviors.

# Parameters
- `index`: the index of the parameter to calculate within `ansatz`
- `ansatz`: the ADAPT state
- `H`: the object defining the cost-function
- `ψ0`: an initial quantum state which the `ansatz` operates on

# Returns
- a number of type `typeof_energy(observable)`.

"""
function ADAPT.partial(
    index::Int,
    ansatz::ADAPT.AbstractAnsatz,
    observable::AnyPauli,
    reference::ADAPT.QuantumState,
)
    state = deepcopy(reference)

    # PARTIAL EVOLUTION
    for i in 1:index
        generator, parameter = ansatz[i]
        ADAPT.evolve_state!(generator, parameter, state)
    end

    # REFLECTION
    generator, _ = ansatz[index]
    costate = -im * generator * state
    # TODO: This line here feels wrong. I'm not sure why.

    # FINISH EVOLUTION
    for i in 1+index:length(ansatz)
        generator, parameter = ansatz[i]
        ADAPT.evolve_state!(generator, parameter, state)
        ADAPT.evolve_state!(generator, parameter, costate)
    end

    return 2 * real(PauliOperators.braket(observable, costate, state))
end


##########################################################################################
#= Default scoring. Variant ADAPT protocols will probably want to override this. =#

function ADAPT.calculate_score(
    ansatz::ADAPT.AbstractAnsatz,
    ::ADAPT.AdaptProtocol,
    generator::AnyPauli,
    observable::AnyPauli,
    reference::ADAPT.QuantumState,
)
    state = ADAPT.evolve_state(ansatz, reference)
    return abs(PauliOperators.measure_commutator(generator, observable, state))
end

function ADAPT.calculate_scores(
    ansatz::ADAPT.AbstractAnsatz,
    ::ADAPT.AdaptProtocol,
    pool::AnyPauli,
    observable::AnyPauli,
    reference::ADAPT.QuantumState,
)
    state = ADAPT.evolve_state(ansatz, reference)
    scores = Vector{typeof_score(ADAPT)}(undef, length(pool))
    for i in eachindex(pool)
        scores[i] = abs(PauliOperators.measure_commutator(pool[i], observable, state))
    end
    return scores
end