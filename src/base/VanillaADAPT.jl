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
    ε = eps(ADAPT.typeof_score(vanilla))
    if all(score -> abs(score) < ε, scores)
        ADAPT.set_converged!(ansatz, true)
        return false
    end

    # MAKE SELECTION
    selected_index = argmax(scores)
    selected_score = scores[selected_index]
    selected_generator = pool[selected_index]
    selected_parameter = zero(ADAPT.typeof_parameter(ansatz))

    # DEFER TO CALLBACKS
    data = ADAPT.Data(
        :scores => scores,
        :selected_index => selected_index,
        :selected_score => selected_score,
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
