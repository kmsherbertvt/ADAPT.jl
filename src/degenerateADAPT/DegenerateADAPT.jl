import ..ADAPT

"""
    DegenerateADAPT

Score pool operators by their initial gradients if they were to be appended to the ansatz.
Equivalently, score pool operators by the expectation value
    of the commutator of the pool operator with the observable.
In the case where the largest scores (gradients) are degenerate between multiple pool operators, choose the 
operator to append to the ansatz randomly.

"""
struct DegenerateADAPT <: ADAPT.AdaptProtocol end
DEG_ADAPT = DegenerateADAPT()

ADAPT.typeof_score(::DegenerateADAPT) = Float64

function ADAPT.adapt!(
    ansatz::ADAPT.AbstractAnsatz,
    trace::ADAPT.Trace,
    adapt_type::DegenerateADAPT,
    pool::ADAPT.GeneratorList,
    observable::ADAPT.Observable,
    reference::ADAPT.QuantumState,
    callbacks::ADAPT.CallbackList,
)
    # CALCULATE SCORES
    scores = ADAPT.calculate_scores(ansatz, adapt_type, pool, observable, reference)

    # CHECK FOR CONVERGENCE
    ε = eps(ADAPT.typeof_score(adapt_type))
    if all(score -> abs(score) < ε, scores)
        ADAPT.set_converged!(ansatz, true)
        return false
    end

    # MAKE SELECTION
    largest_score = maximum(scores); indices_of_largest_scores = findall(a->a==largest_score, scores) 
#     println("gradient degeneracy: ",length(indices_of_largest_scores))
#     println("operators with degenerate max. gradients: ",pool[indices_of_largest_scores])
    selected_index = rand(indices_of_largest_scores)    
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
        stop = stop || callback(data, ansatz, trace, adapt_type, pool, observable, reference)
        # Note that, as soon as `stop` is true, subsequent callbacks are short-circuited.
    end
    (stop || ADAPT.is_converged(ansatz)) && return false

    # PERFORM ADAPTATION
    push!(ansatz, selected_generator => selected_parameter)
    ADAPT.set_optimized!(ansatz, false)
    return true
end

function ADAPT.calculate_score(
    ansatz::ADAPT.AbstractAnsatz,
    adapt_type::DegenerateADAPT,
    generator::ADAPT.Generator,
    observable::ADAPT.Observable,
    reference::ADAPT.QuantumState,
)
    L = length(ansatz)
    candidate = deepcopy(ansatz)
    push!(candidate, generator => zero(ADAPT.typeof_parameter(ansatz)))
    return abs(ADAPT.partial(L+1, candidate, observable, reference))
end