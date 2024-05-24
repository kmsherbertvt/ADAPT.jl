import ..ADAPT
import PauliOperators: ScaledPauliVector

"""
    TETRISADAPT

Score pool operators by their initial gradients if they were to be appended to the pool.
TETRIS-ADAPT is a modified version of ADAPT-VQE in which multiple operators with disjoint 
support are added to the ansatz at each iteration. They are chosen by selecting from 
operators ordered in decreasing magnitude of gradients.
"""
struct TETRISADAPT <: ADAPT.AdaptProtocol end
TETRIS = TETRISADAPT()

ADAPT.typeof_score(::TETRISADAPT) = Float64

function support(spv::ScaledPauliVector)
    indices = Vector{Int64}()
    for sp in spv
        op_indices=findall(x -> x != 'I',string(sp.pauli))
        append!(indices,op_indices)
    end
    indices = unique(indices)
    return indices
end

function ADAPT.adapt!(
    ansatz::ADAPT.AbstractAnsatz,
    trace::ADAPT.Trace,
    adapt_type::TETRISADAPT,
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
    #= TODO: there is probably a better way to make the selection of operators, given possible 
    gradient degeneracies =#
    # for now, we start by choosing the op. with the largest gradient (in the case of 
    # degeneracy, argmax gives the first index)
    scores_copy = deepcopy(scores);
    largest_score = maximum(scores_copy); index_of_largest_score = argmax(scores_copy) 
    ops_to_add = [index_of_largest_score]
    current_support = support(pool[ops_to_add[1]]); max_support = length(string(pool[1][1].pauli))
    for count in range(1,length(scores_copy)-1)
        for i in ops_to_add
            scores_copy[i] = 0.0
        end
        next_largest_score = maximum(scores_copy); index_of_tentative_op = argmax(scores_copy) 
        if next_largest_score > 1e-3
            op_is_disjoint = true
            if !isdisjoint(support(pool[index_of_tentative_op]), current_support)
                op_is_disjoint = false
            end
            if op_is_disjoint
                push!(ops_to_add, index_of_tentative_op)
                append!(current_support, support(pool[index_of_tentative_op]))
            end
        end
        scores_copy[index_of_tentative_op] = 0.0
        if length(current_support) == max_support || next_largest_score < 1e-3
            break
        end
    end  
    selected_indices = ops_to_add
    selected_scores = scores[selected_indices];
    selected_generators = pool[selected_indices];
    selected_parameters = zeros(ADAPT.typeof_parameter(ansatz), length(ops_to_add));

    # DEFER TO CALLBACKS
    data = ADAPT.Data(
        :scores => scores,
        :selected_index => selected_indices,
        :selected_score => selected_scores,
        :selected_generator => selected_generators,
        :selected_parameter => selected_parameters,
    )

    stop = false
    for callback in callbacks
        stop = stop || callback(data, ansatz, trace, adapt_type, pool, observable, reference)
        # Note that, as soon as `stop` is true, subsequent callbacks are short-circuited.
    end
    (stop || ADAPT.is_converged(ansatz)) && return false

    # PERFORM ADAPTATION
    for i in range(1,length(selected_generators))
        push!(ansatz, selected_generators[i] => selected_parameters[i])
    end
    ADAPT.set_optimized!(ansatz, false)
    return true
end

function ADAPT.calculate_score(
    ansatz::ADAPT.AbstractAnsatz,
    adapt_type::TETRISADAPT,
    generator::ADAPT.Generator,
    observable::ADAPT.Observable,
    reference::ADAPT.QuantumState,
)
    L = length(ansatz)
    candidate = deepcopy(ansatz)
    push!(candidate, generator => zero(ADAPT.typeof_parameter(ansatz)))
    return abs(ADAPT.partial(L+1, candidate, observable, reference))
end