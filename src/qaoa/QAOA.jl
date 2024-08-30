import ..ADAPT

"""
    QAOAAnsatz{F<:Parameter,G<:Generator}(
        observable::G,
        γ0::F,
        generators::Vector{G},
        β_parameters::Vector{F},
        γ_parameters::Vector{F},
        optimized::Bool,
        converged::Bool,
    )

An ADAPT state suitable for ADAPT-QAOA.
The standard ADAPT generators are interspersed with the observable itself.

# Type Parameters
- `F`: the number type for the parameters (usually `Float64` is appropriate).
- `G`: the generator type. Uniquely for QAOA, `G` must ALSO be a valid `Observable` type.

# Parameter
- `observable`: the observable, which is interspersed with generators when evolving
- `γ0`: initial coefficient of the observable, whenever a new generator is added
- `generators`: list of current generators (i.e. mixers)
- `β_parameters`: list of current generator coefficients
- `γ_parameters`: list of current observable coefficients
- `optimized`: whether the current parameters are flagged as optimal
- `converged`: whether the current generators are flagged as converged

"""
struct QAOAAnsatz{F,G} <: ADAPT.AbstractAnsatz{F,G}
    observable::G
    γ0::F
    generators::Vector{G}
    β_parameters::Vector{F}
    γ_parameters::Vector{F}
    optimized::Ref{Bool}
    converged::Ref{Bool}
end

ADAPT.__get__generators(ansatz::QAOAAnsatz) = vec(permutedims(
    hcat(fill(ansatz.observable, length(ansatz.generators)), ansatz.generators)
))
ADAPT.__get__parameters(ansatz::QAOAAnsatz) = vec(permutedims(
    hcat(ansatz.γ_parameters, ansatz.β_parameters)
))
    #= TODO: These functions are actually redundant.
        We should probably alter the `AbstractAnsatz` interface. =#

ADAPT.__get__optimized(ansatz::QAOAAnsatz) = ansatz.optimized
ADAPT.__get__converged(ansatz::QAOAAnsatz) = ansatz.converged

"""
    Ansatz(γ0, observable)

Convenience constructor for initializing an empty ansatz.

# Parameters
- γ0
- observable

Note that, uniquely for QAOA,
    the observable and the pool operators must be of the same type.

"""
QAOAAnsatz(γ0, observable) = QAOAAnsatz(
    observable, γ0,
    typeof(observable)[],
    typeof(γ0)[],
    typeof(γ0)[],
    Ref(true), Ref(false),
)


##########################################################################################
#= AbstractVector interface. =#

#= Convenience function for extracting the half-index of integer `ix`. =#
half(ix) = 1 + ((ix-1) >> 1)

Base.size(ansatz::QAOAAnsatz) = size(ansatz.generators) .<< 1
Base.IndexStyle(::Type{<:QAOAAnsatz}) = IndexLinear()

function Base.getindex(ansatz::QAOAAnsatz, i::Int)
    ((i-1) & 1 == 0) && return ansatz.observable => ansatz.γ_parameters[half(i)]
    return ansatz.generators[half(i)] => ansatz.β_parameters[half(i)]
end

function Base.setindex!(ansatz::QAOAAnsatz{F,G}, pair::Pair{G,F}, i::Int) where {F,G}
    #= TODO: We are making a major assumption,
            that setindex! is only called in the context of push!(G => x),
            i.e. attaching a new generator.
        Thus, we assume `ansatz[i] = (G => x)` is never called. =#
    ansatz.generators[half(i)] = pair.first
    ansatz.β_parameters[half(i)] = pair.second
    ansatz.γ_parameters[half(i)] = ansatz.γ0
end

function Base.resize!(ansatz::QAOAAnsatz, nl::Int)
    resize!(ansatz.generators, half(nl))
    resize!(ansatz.β_parameters, half(nl))
    resize!(ansatz.γ_parameters, half(nl))
end

##########################################################################################
#= Evolution. =#

function ADAPT.angles(ansatz::QAOAAnsatz)
    #= TODO: We should be able to make this a view, to avoid allocations.

    Seems we ought replace a whole bunch of `copy(angles...)` with `collect(angles...)`,
        but otherwise such an implementation should be fine.

    =#
    return vec(permutedims(hcat(ansatz.γ_parameters, ansatz.β_parameters)))
end

function ADAPT.bind!(ansatz::QAOAAnsatz{F,G}, x::AbstractVector{F}) where {F,G}
    x = reshape(x, 2, :)
    ansatz.γ_parameters .= @view(x[1,:])
    ansatz.β_parameters .= @view(x[2,:])
end

##########################################################################################
#= Improved scoring for vanilla ADAPT. =#

#= TODO: There's a major dispatch ambiguity problem.

For now I'm hard-coding the generator/observable/protocol types,
    but this really shouldn't be necessary.

I think adding a `calculate_score(..., state)` into the core interface,
    so that `calculate_score(ansatz, ..., reference)` has an obvious dispatch
    (and so does `calculate_scores(...)`, for that matter!)
    may help, but I'm not sure it solves it.

Need to think it through a bit more.

But ultimately we shouldn't need to import PauliOperators here.

=#

import ADAPT.Basics.MyPauliOperators
import ADAPT.Basics.MyPauliOperators: Pauli, ScaledPauli, PauliSum, ScaledPauliVector
AnyPauli = Union{Pauli, ScaledPauli, PauliSum, ScaledPauliVector}
# TODO: Replace `MyPauliOperators` with `PauliOperators` throughout, once merged.

function ADAPT.calculate_score(
    ansatz::QAOAAnsatz,
    ::ADAPT.Basics.VanillaADAPT,
    generator::AnyPauli,
    observable::AnyPauli,
    reference::ADAPT.QuantumState,
)
    state = ADAPT.evolve_state(ansatz, reference)
    ADAPT.evolve_state!(ansatz.observable, ansatz.γ0, state)
    return abs(ADAPT.Basics.MyPauliOperators.measure_commutator(
            generator, observable, state))
end

function ADAPT.calculate_score(
    ansatz::QAOAAnsatz,
    ::ADAPT.Degenerate_ADAPT.DegenerateADAPT,
    generator::AnyPauli,
    observable::AnyPauli,
    reference::ADAPT.QuantumState,
)
    state = ADAPT.evolve_state(ansatz, reference)
    ADAPT.evolve_state!(ansatz.observable, ansatz.γ0, state)
    return abs(ADAPT.Basics.MyPauliOperators.measure_commutator(
            generator, observable, state))
end

##########################################################################################
#= HACK - Assume that ScaledPauliVector generators appearing in QAOA
            consist only of commuting terms.

    Motivation: evolve the Hamiltonian (explicitly commuting diagonal terms) faster.

=#

import LinearAlgebra: dot

"""

Carbon copy of the usual costate function with Pauli operators.

The only difference is that we don't copy in the method defining special behavior
    for ScaledPauliVectors in this namespace,
    so those will be treated as though they consist only of commuting terms.

"""
function __make__costate(G, x, Ψ)
    costate = -im * G * Ψ
    ADAPT.evolve_state!(G, x, costate)
    return costate
end

"""

Carbon copy of the usual gradient with Pauli operators.

The only difference is which `__make__costate` function is getting called.

"""
function ADAPT.gradient!(
    result::AbstractVector,
    ansatz::QAOAAnsatz,
    observable::AnyPauli,
    reference::ADAPT.QuantumState,
)
    ψ = ADAPT.evolve_state(ansatz, reference)   # FULLY EVOLVED ANSATZ |ψ⟩
    λ = observable * ψ                          # CALCULATE |λ⟩ = H |ψ⟩

    for i in reverse(eachindex(ansatz))
        G, θ = ansatz[i]
        ADAPT.evolve_state!(G', -θ, ψ)          # UNEVOLVE BRA
        σ = __make__costate(G, θ, ψ)            # CALCULATE |σ⟩ = exp(-iθG) (-iG) |ψ⟩
        result[i] = 2 * real(dot(σ, λ))         # CALCULATE GRADIENT ⟨λ|σ⟩ + h.t.
        ADAPT.evolve_state!(G', -θ, λ)          # UNEVOLVE KET
    end

    return result
end

# TODO: In principle we could get even faster by replacing the `evolve_state!` above with a local version that assumes diagonality. This would save a lot of memory allocation so I anticipate a noticable time improvement, but it would not be asymptotically faster, so I don't deem it worth the effort. (If we do do this, also override evolve_state!(ansatz, state) and call the same method.)