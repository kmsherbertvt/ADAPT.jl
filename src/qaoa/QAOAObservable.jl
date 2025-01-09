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
# TODO: Major dispatch ambiguity problem, see a little more at the bottom of the file.
ADAPT.calculate_score(
    ansatz::ADAPT.AbstractAnsatz,
    adapt::ADAPT.Degenerate_ADAPT.DegenerateADAPT,
    generator::ADAPT.Basics.AnyPauli,
    dspv::QAOAObservable,
    reference::ADAPT.QuantumState,
) = ADAPT.calculate_score(ansatz, adapt, generator, dspv.spv, reference)
ADAPT.calculate_scores(
    ansatz::ADAPT.AbstractAnsatz,
    adapt::ADAPT.Degenerate_ADAPT.DegenerateADAPT,
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
