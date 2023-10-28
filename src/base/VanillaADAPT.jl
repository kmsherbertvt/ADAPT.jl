import ..ADAPT

"""
    VanillaADAPT

Score pool operators by their initial gradients if they were to be appended to the pool.
Equivalently, score pool operators by the expectation value
    of the commutator of the pool operator with the observable.

This protocol is defined for when the pool operators and the observable are AbstractPaulis.
Note that fermionic operators are perfectly well-represented with AbstractPaulis.

"""
struct VanillaADAPT <: ADAPT.AdaptProtocol end
VANILLA = VanillaADAPT()

ADAPT.typeof_score(::VanillaADAPT) = Float64
#= NOTE: Strict typing is consequent of strict typing in the PauliOperators package. =#

"""
    measure_commutator(
        A::AbstractPauli,
        B::AbstractPauli,
        Ψ::Union{SparseKetBasis,AbstractVector},
    )

Calculate the expectation value of the commutator, ie. ⟨Ψ|[A,B]|Ψ⟩.

TODO: There *could* be a place for this in PauliOperators,
        but it would need to be carefully fleshed out type by type.
    A and B needn't be Hermitian in general (though I assume they are here),
        so my intuition is rather lacking.

"""
function measure_commutator(
    A::AbstractPauli,
    B::AbstractPauli,
    Ψ::Union{SparseKetBasis,AbstractVector},
)
    commutator = A*B - B*A
    #= TODO: ALLOCATIONS!
        Diksha has a much more efficient version in the ACSE package.
        But, it would need distinct methods for all the type combinations.
        All my operators are fixed size,
            so I don't know that I care about the allocations here.
    =#
    return expectation_value(commutator, Ψ)
end

function ADAPT.calculate_score(
    ansatz::ADAPT.AbstractAnsatz,
    ::VanillaADAPT,
    generator::AbstractPauli,
    observable::AbstractPauli,
    reference::ADAPT.QuantumState,
)
    state = ADAPT.evolve_state(ansatz, reference)
    return real(measure_commutator(generator, observable, state))
end

function ADAPT.calculate_scores(
    ansatz::ADAPT.AbstractAnsatz,
    ::VanillaADAPT,
    pool::AbstractPauli,
    observable::AbstractPauli,
    reference::ADAPT.QuantumState,
)
    state = ADAPT.evolve_state(ansatz, reference)
    scores = Vector{typeof_score(ADAPT)}(undef, length(pool))
    for i in eachindex(pool)
        scores[i] = real(measure_commutator(pool[i], observable, state))
    end
    return scores
end

function ADAPT.adapt!(
    ansatz::ADAPT.AbstractAnsatz,
    trace::ADAPT.Trace,
    vanilla::VanillaADAPT,
    pool::ADAPT.GeneratorList,
    observable::ADAPT.Observable,
    reference::ADAPT.QuantumState,
    callbacks::ADAPT.CallbackList,
)
    # CALCULATE SCORES
    scores = ADAPT.calculate_scores(ansatz, vanilla, pool, observable, reference)

    # CHECK FOR CONVERGENCE
    ε = eps(typeof_score(vanilla))
    if any(score -> abs(score) ≥ ε, scores)
        ADAPT.set_converged!(ansatz, true)
        return false
    end

    # MAKE SELECTION
    selected_index = argmax(scores)
    selected_generator = pool[pool_index]
    selected_parameter = zero(typeof_parameter(ansatz))

    # DEFER TO CALLBACKS
    data = Data(
        :scores => scores,
        :selected_index => selected_index,
        :selected_generator => selected_generator,
        :selected_parameter => selected_parameter,
    )

    stop = false
    for callback in callbacks
        stop = stop || callback(data, ansatz, trace, vanilla, pool, observable, reference)
        # Note that, as soon as `stop` is true, subsequent callbacks are short-circuited.
    end
    (stop || ADAPT.is_converged(ansatz)) && return false

    # PERFORM ADAPTATION
    push!(ansatz, selected_generator => selected_parameter)
    ADAPT.set_optimized!(ansatz, false)
    return true
end
