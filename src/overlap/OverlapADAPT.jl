import ..ADAPT
import PauliOperators: SparseKetBasis

import LinearAlgebra: dot

struct Infidelity{T<:ADAPT.QuantumState}
    Φ::T
end

"""
    Matrix(infidelity)

Convert an infidelity to a matrix.

This implementation assumes:
- The target state `infidelity.Φ` can be cast to a vector.
- The reference state in `evaluate(infidelity, reference)` is always normalized.

"""
function Base.Matrix(infidelity::Infidelity)
    Φ = Vector(infidelity.Φ)
    ρ = Φ*Φ'
    return one(ρ) - ρ
end

function ADAPT.typeof_energy(observable::Infidelity)
    return real(eltype(observable.Φ))
end

function ADAPT.typeof_energy(observable::Infidelity{SparseKetBasis{N,T}}) where {N,T}
    return real(T)
end

function ADAPT.evaluate(observable::Infidelity, Ψ::ADAPT.QuantumState)
    return 1 - abs2(dot(observable.Φ, Ψ))
end

function ADAPT.partial(
    index::Int,
    ansatz::ADAPT.AbstractAnsatz,
    observable::Infidelity,
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
    costate = ADAPT.Basics.__make__costate(generator, parameter, state)
    ADAPT.evolve_state!(generator, parameter, state)

    # FINISH EVOLUTION
    for i in 1+index:length(ansatz)
        generator, parameter = ansatz[i]
        ADAPT.evolve_state!(generator, parameter, state)
        ADAPT.evolve_state!(generator, parameter, costate)
    end

    return -2 * real(dot(state, observable.Φ) * dot(observable.Φ, costate))
end

function ADAPT.gradient!(
    result::AbstractVector,
    ansatz::ADAPT.AbstractAnsatz,
    observable::Infidelity,
    reference::ADAPT.QuantumState,
)
    ψ = ADAPT.evolve_state(ansatz, reference)           # FULLY EVOLVED ANSATZ |ψ⟩
    λ = deepcopy(observable.Φ) * dot(observable.Φ, ψ)   # CALCULATE |λ⟩ = |ϕ⟩⟨ϕ|ψ⟩

    for i in reverse(eachindex(ansatz))
        G, θ = ansatz[i]
        ADAPT.evolve_state!(G', -θ, ψ)          # UNEVOLVE BRA
        σ = ADAPT.Basics.__make__costate(G, θ, ψ)   # CALCULATE |σ⟩ = exp(-iθG) (-iG) |ψ⟩
        result[i] = -2 * real(dot(σ, λ))        # CALCULATE GRADIENT -(⟨λ|σ⟩ + h.t.)
        ADAPT.evolve_state!(G', -θ, λ)          # UNEVOLVE KET
    end

    return result
end

function ADAPT.calculate_score(
    ansatz::ADAPT.AbstractAnsatz,
    ::ADAPT.AdaptProtocol,
    generator::ADAPT.Generator,
    observable::Infidelity,
    reference::ADAPT.QuantumState,
)
    state = ADAPT.evolve_state(ansatz, reference)
    x = dot(observable.Φ, state)
    x̄ = dot(generator * observable.Φ, state)
    return 2 * abs(real(x̄)*imag(x) - imag(x̄)*real(x))
end

function ADAPT.calculate_scores(
    ansatz::ADAPT.AbstractAnsatz,
    adapt::ADAPT.AdaptProtocol,
    pool::ADAPT.GeneratorList,
    observable::Infidelity,
    reference::ADAPT.QuantumState,
)
    state = ADAPT.evolve_state(ansatz, reference)
    x = dot(observable.Φ, state)
    scores = Vector{ADAPT.typeof_score(adapt)}(undef, length(pool))
    for i in eachindex(pool)
        x̄ = dot(pool[i] * observable.Φ, state)
        scores[i] = 2 * abs(real(x̄)*imag(x) - imag(x̄)*real(x))
    end
    return scores
end