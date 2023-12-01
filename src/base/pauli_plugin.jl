import ..ADAPT

import .MyPauliOperators
import .MyPauliOperators: SparseKetBasis
import .MyPauliOperators: FixedPhasePauli, ScaledPauli, Pauli
import .MyPauliOperators: ScaledPauliVector, PauliSum
PauliOperators = MyPauliOperators
# TODO: Replace `MyPauliOperators` with `PauliOperators`

import LinearAlgebra: mul!

AnyPauli = Union{Pauli, ScaledPauli, PauliSum, ScaledPauliVector}
#= NOTE: ANY Pauli is a little strong.
We don't actually support FixedPhasePauli in this library.

(This is just because it's rather confusingly defined,
    so I want users to be explicit and just use Pauli.
    I bet that's what Nick intended also.)
=#

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
    G::Pauli,
    θ::ADAPT.Parameter,
    Ψ::ADAPT.QuantumState,
)
    angle = -θ * PauliOperators.get_phase(G) * PauliOperators.get_phase(G.pauli)'
    PauliOperators.cis!(Ψ, G.pauli, angle)
    Ψ isa SparseKetBasis && PauliOperators.clip!(Ψ)
end

function ADAPT.evolve_state!(
    G::ScaledPauli,
    θ::ADAPT.Parameter,
    Ψ::ADAPT.QuantumState,
)
    angle = -θ * G.coeff * PauliOperators.get_phase(G.pauli)'
    PauliOperators.cis!(Ψ, G.pauli, angle)
    Ψ isa SparseKetBasis && PauliOperators.clip!(Ψ)
end

function ADAPT.evolve_state!(
    G::ScaledPauliVector,
    θ::ADAPT.Parameter,
    Ψ::ADAPT.QuantumState,
)
    for P in G
        ADAPT.evolve_state!(P, θ, Ψ)
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
    for i in 1:index-1
        generator, parameter = ansatz[i]
        ADAPT.evolve_state!(generator, parameter, state)
    end

    # REFLECTION
    generator, parameter = ansatz[index]
    costate = __make__costate(generator, parameter, state)
    ADAPT.evolve_state!(generator, parameter, state)

    # FINISH EVOLUTION
    for i in 1+index:length(ansatz)
        generator, parameter = ansatz[i]
        ADAPT.evolve_state!(generator, parameter, state)
        ADAPT.evolve_state!(generator, parameter, costate)
    end

    return 2 * real(PauliOperators.braket(observable, costate, state))
end


"""
    __make__costate(G, x, Ψ)

Compute ∂/∂x exp(ixG) |ψ⟩.

"""
function __make__costate(G, x, Ψ)
    costate = -im * G * Ψ
    ADAPT.evolve_state!(G, x, costate)
    return costate
end

"""
    __make__costate(G::ScaledPauliVector, x, Ψ)

Compute ∂/∂x exp(ixG) |ψ⟩.

Default implementation just applies -iG to Ψ then evolves.
That's fine as long as the evolution is exact.
But evolution is not exact if `G` is a `ScaledPauliVector` containing non-commuting terms.
In such a case, the co-state must be more complicated.

"""
function __make__costate(G::ScaledPauliVector, x, Ψ::SparseKetBasis)
    costate = zero(Ψ)
    for (i, P) in enumerate(G)
        term = deepcopy(Ψ)
        for j in 1:i-1; ADAPT.evolve_state!(G[j], x, term); end         # RIGHT EVOLUTION
        term = -im * P * term                                           # REFLECTION
        #= TODO: Hypothetically could implement mul! for SparseKetBasis someday? =#
        for j in i:length(G); ADAPT.evolve_state!(G[j], x, term); end   # LEFT EVOLUTION
        sum!(costate, term)
    end
    PauliOperators.clip!(costate)
    return costate
end

function __make__costate(G::ScaledPauliVector, x, Ψ::AbstractVector)
    costate = zero(Ψ)
    work_l = Array{eltype(Ψ)}(undef, size(Ψ))
    work_r = Array{eltype(Ψ)}(undef, size(Ψ))
    for (i, P) in enumerate(G)
        work_r .= Ψ
        for j in 1:i-1; ADAPT.evolve_state!(G[j], x, work_r); end       # RIGHT EVOLUTION
        mul!(work_l, P, work_r); work_l .*= -im                         # REFLECTION
        for j in i:length(G); ADAPT.evolve_state!(G[j], x, work_l); end # LEFT EVOLUTION
        costate .+= work_l
    end
    return costate
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