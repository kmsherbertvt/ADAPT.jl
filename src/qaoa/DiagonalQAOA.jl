import ..ADAPT
import PauliOperators: ScaledPauliVector, SparseKetBasis

"""
    QAOAObservable(spv::ScaledPauliVector)

Wrap a ScaledPauliVector observable in a view that assumes each element is diagonal,
    allowing for more memory-efficient state evolution.

The constructor throws an error if any element `sp` of `spv` has `sp.pauli.x != 0`.

"""
struct QAOAObservable{N}
    spv::ScaledPauliVector{N}

    function QAOAObservable(spv::ScaledPauliVector{N}) where {N}
        for sp in spv
            @assert sp.pauli.x == 0
        end
        return new{N}(spv)
    end
end



function ADAPT.evolve_state!(
    dspv::QAOAObservable,
    θ::ADAPT.Parameter,
    Ψ::SparseKetBasis,
)
    for ket in keys(Ψ)
        α = zero(θ)                 # TOTAL PHASE SHIFT FOR THIS KET
        for sp in dspv.spv
            phase = iseven(count_ones(sp.pauli.z & ket.v)) ? 1 : -1
            α += phase * sp.coeff
        end
        Ψ[ket] *= exp(-im * α * θ)
    end
    Ψ
end

function ADAPT.evolve_state!(
    dspv::QAOAObservable,
    θ::ADAPT.Parameter,
    Ψ::AbstractVector,
)
    for i in eachindex(Ψ)
        v = i - 1
        α = zero(θ)                 # TOTAL PHASE SHIFT FOR THIS KET
        for sp in dspv.spv
            phase = iseven(count_ones(sp.pauli.z & v)) ? 1 : -1
            α += phase * sp.coeff
        end
        Ψ[i] *= exp(-im * α * θ)
    end
    return Ψ
end


# Delegate to ScaledPauliVector for unhandled state types.
ADAPT.evolve_state!(
    dspv::QAOAObservable,
    θ::ADAPT.Parameter,
    Ψ::ADAPT.QuantumState,
) = ADAPT.evolve_state!(dspv.spv, θ, Ψ)

# Delegate all the Observable behavior directly to the ScaledPauliVector.
ADAPT.evaluate(
    dspv::QAOAObservable,
    state::ADAPT.QuantumState,
) = ADAPT.evaluate(dspv.spv, state)
ADAPT.partial(
    index::Int,
    ansatz::ADAPT.AbstractAnsatz,
    dspv::QAOAObservable,
    reference::ADAPT.QuantumState,
) = ADAPT.partial!(index, ansatz, dspv.spv, reference)
ADAPT.gradient!(
    result::AbstractVector,
    ansatz::ADAPT.AbstractAnsatz,
    dspv::QAOAObservable,
    reference::ADAPT.QuantumState,
) = ADAPT.gradient!(result, ansatz, dspv.spv, reference)
ADAPT.calculate_score(
    ansatz::ADAPT.AbstractAnsatz,
    adapt::ADAPT.Basics.VanillaADAPT,
    generator::ADAPT.Basics.AnyPauli,
    dspv::QAOAObservable,
    reference::ADAPT.QuantumState,
) = ADAPT.calculate_score(ansatz, adapt, generator, dspv.spv, reference)
ADAPT.calculate_scores(
    ansatz::ADAPT.AbstractAnsatz,
    adapt::ADAPT.Basics.VanillaADAPT,
    pool::ADAPT.Basics.AnyPauli,
    dspv::QAOAObservable,
    reference::ADAPT.QuantumState,
) = ADAPT.calculate_scores(ansatz, adapt, pool, dspv.spv, reference)

# Delegate multiplication (for co-state construction in gradient)
Base.:*(dspv::QAOAObservable, factor) = dspv.spv * factor
Base.:*(factor, dspv::QAOAObservable) = factor * dspv.spv
Base.:*(
    dspv1::QAOAObservable,
    dspv2::QAOAObservable,
) = QAOAObservable(dspv1.spv * dspv2.spv)

# Delegate adjoint (for reverse evolution in gradient)
Base.adjoint(dspv::QAOAObservable) = QAOAObservable(adjoint(dspv.spv))

# Delegate matrix construction (for convenient exact diagonalization).
Base.Matrix(dspv::QAOAObservable) = Matrix(dspv.spv)

#= NOTE: The most important thing to the gradient calculation
    is using the version of `Basics.make_costate` that assumes all terms commute
    (MUCH more efficient).
The way the `pauli_plugin` is currently written, this is the default,
    so it doesn't need to be made explicit here. =#



"""
    DiagonalQAOAAnsatz{F<:Parameter,G<:Generator}(
        observable::QAOAObservable,
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
- `G`: the generator type.

# Parameter
- `observable`: the observable, which is interspersed with generators when evolving
- `γ0`: initial coefficient of the observable, whenever a new generator is added
- `generators`: list of current generators (i.e. mixers)
- `β_parameters`: list of current generator coefficients
- `γ_parameters`: list of current observable coefficients
- `optimized`: whether the current parameters are flagged as optimal
- `converged`: whether the current generators are flagged as converged

"""
struct DiagonalQAOAAnsatz{F,G} <: ADAPT.AbstractAnsatz{F,G}
    observable::QAOAObservable
    γ0::F
    generators::Vector{G}
    β_parameters::Vector{F}
    γ_parameters::Vector{F}
    optimized::Ref{Bool}
    converged::Ref{Bool}
end

ADAPT.__get__generators(ansatz::DiagonalQAOAAnsatz) = vec(permutedims(
    hcat(fill(ansatz.observable, length(ansatz.generators)), ansatz.generators)
))
ADAPT.__get__parameters(ansatz::DiagonalQAOAAnsatz) = vec(permutedims(
    hcat(ansatz.γ_parameters, ansatz.β_parameters)
))
    #= TODO: These functions are actually redundant.
        We should probably alter the `AbstractAnsatz` interface. =#

ADAPT.__get__optimized(ansatz::DiagonalQAOAAnsatz) = ansatz.optimized
ADAPT.__get__converged(ansatz::DiagonalQAOAAnsatz) = ansatz.converged

"""
    Ansatz(γ0, observable)

Convenience constructor for initializing an empty ansatz.

# Parameters
- γ0
- observable

Note that, uniquely for QAOA,
    the observable and the pool operators must be of the same type.

"""
DiagonalQAOAAnsatz(γ0, pool, observable) = DiagonalQAOAAnsatz(
    observable, γ0,
    eltype(pool)[],
    typeof(γ0)[],
    typeof(γ0)[],
    Ref(true), Ref(false),
)


##########################################################################################
#= AbstractVector interface. =#

#= Convenience function for extracting the half-index of integer `ix`. =#
# half(ix) = 1 + ((ix-1) >> 1)

Base.size(ansatz::DiagonalQAOAAnsatz) = size(ansatz.generators) .<< 1
Base.IndexStyle(::Type{<:DiagonalQAOAAnsatz}) = IndexLinear()

function Base.getindex(ansatz::DiagonalQAOAAnsatz, i::Int)
    ((i-1) & 1 == 0) && return ansatz.observable => ansatz.γ_parameters[half(i)]
    return ansatz.generators[half(i)] => ansatz.β_parameters[half(i)]
end

function Base.setindex!(ansatz::DiagonalQAOAAnsatz{F,G}, pair::Pair{G,F}, i::Int) where {F,G}
    #= TODO: We are making a major assumption,
            that setindex! is only called in the context of push!(G => x),
            i.e. attaching a new generator.
        Thus, we assume `ansatz[i] = (G => x)` is never called. =#
    ansatz.generators[half(i)] = pair.first
    ansatz.β_parameters[half(i)] = pair.second
    ansatz.γ_parameters[half(i)] = ansatz.γ0
end

function Base.resize!(ansatz::DiagonalQAOAAnsatz, nl::Int)
    resize!(ansatz.generators, half(nl))
    resize!(ansatz.β_parameters, half(nl))
    resize!(ansatz.γ_parameters, half(nl))
end

##########################################################################################
#= Evolution. =#

function ADAPT.angles(ansatz::DiagonalQAOAAnsatz)
    #= TODO: We should be able to make this a view, to avoid allocations.

    Seems we ought replace a whole bunch of `copy(angles...)` with `collect(angles...)`,
        but otherwise such an implementation should be fine.

    =#
    return vec(permutedims(hcat(ansatz.γ_parameters, ansatz.β_parameters)))
end

function ADAPT.bind!(ansatz::DiagonalQAOAAnsatz{F,G}, x::AbstractVector{F}) where {F,G}
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
    ansatz::DiagonalQAOAAnsatz,
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
    ansatz::DiagonalQAOAAnsatz,
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