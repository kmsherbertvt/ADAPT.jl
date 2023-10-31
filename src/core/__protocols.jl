#= DISCLAIMER: This file is meant to be included in ADAPT.jl.
    Names not defined in this file or in base Julia can be found there.
=#

##########################################################################################
#= Data passing interface. =#

"""
    Trace

Semantic alias for a compact record of the entire ADAPT run.

The `adapt!`, `optimize!`, and `run!` functions require a `Trace` object,
    which will be mutated throughout according to callbacks.
Initialize an empty trace object by `trace = Trace()`.

Keys can in principle be any Symbol at all.
You can design your own protocols to fill the data,
    and your own callbacks to use it.
That said, see the `Callbacks` module for some standard choices.

"""
Trace = Dict{Symbol, Any}

"""
    Data

Semantic alias for trace-worthy information from a single adapt or vqe iteration.

You'll never actually have to deal with this object
    unless you are implementing your own protocol or callback.

Keys can in principle be any Symbol at all.
You can design your own protocols to fill the data,
    and your own callbacks to use it.
That said, see the `Callbacks` module for some standard choices.

"""
Data = Dict{Symbol, Any}

"""
    AbstractCallback

A function to be called at each adapt iteration, or vqe iteration, or both.

# Common Examples
1. Tracers: update the running `trace` with information passed in `data`
2. Printers: display the information passed in `data` to the screen or to a file
3. Stoppers: flag the ADAPT state as converged, based on some condition

In particular, the standard way to converge an adapt run is to include a `ScoreStopper`.
Otherwise, the run will keep adding parameters until every score is essentially zero.

More details can be found in the `Callbacks` module,
    where many standard callbacks are already implemented.

# Implementation

Callbacks are implemented as callable objects, with two choices of method header
    (one for adaptations, one for optimization iterations).

- `(::AbstractCallback)(
    ::Data,
    ::AbstractAnsatz,
    ::Trace,
    ::AdaptProtocol,
    ::GeneratorList,
    ::Observable,
    ::QuantumState,
  )`

- `(::AbstractCallback)(
    ::Data,
    ::AbstractAnsatz,
    ::Trace,
    ::OptimizationProtocol,
    ::Observable,
    ::QuantumState,
  )`

If your callback is only meant for adaptations,
    simply do not implement the method for optimizations.
(Behind the scenes, every `AbstractCallback` has default implementations for both methods,
    which just don't do anything.)

Precisely what data is contained within the `data` depends on the protocol.
For example, the `ScoreStopper` expects to find the key `:scores`,
    whose value is a `ScoreList`, one score for each pool operator.
Generally, the callback should assume `data` has whatever it needs,
    and if it doesn't, that means this callback is incompatible with the given protocol.
That said, see the `Callbacks` module for some standard choices.

The callback is free to mutate the `ansatz`.
For example, the `ScoreStopper` signals a run should end by calling `set_converged!`.
But, if the callback wants to signal the run should end in an UN-converged state,
    it should simply return `true`.

"""
abstract type AbstractCallback end

"""
    CallbackList

Semantic alias for a vector of callbacks.

"""
CallbackList = AbstractVector{<:AbstractCallback}


##########################################################################################
#= Adapt protocol. =#

"""
    AdaptProtocol

A distinctive protocol for adding new parameters after optimizing an initial ansatz.

# Implementation

Sub-types must implement the following method:
- `typeof_parameter(::AbstractAnsatz)::Type{<:Parameter}`

In addition, there must be a compatible implementation for:
- `calculate_score(
    ::AbstractAnsatz,
    ::AdaptProtocol,
    ::Generator,
    ::Observable,
    ::QuantumState,
  )::Score`

- `adapt!(
    ::AbstractAnsatz,
    ::Trace,
    ::AdaptProtocol,
    ::GeneratorList,
    ::Observable,
    ::QuantumState,
    ::CallbackList,
  )`

For the most part, sub-types should be singleton objects, ie. no attributes.
Arbitrary hyperparameters like gradient tolerance should be delegated to callbacks
    as much as possible.
It's okay to insist a particular callback always be included with your ADAPT protocol,
    so long as you are clear in the documentation.
That said, this rule is a "style" guideline, not a contract.

"""
abstract type AdaptProtocol end


"""
    typeof_score(::AdaptProtocol)

The number type of the score for each pool operator.

# Implementation

Usually a sub-type of `AbstractFloat`, and probably just about always `Float64`.

"""
typeof_score(::AdaptProtocol) = NotImplementedError()

"""
    calculate_score(
        ansatz::AbstractAnsatz,
        ADAPT::AdaptProtocol,
        generator::Generator,
        H::Observable,
        ψ0::QuantumState,
    )

Calculate an "importance" score for a generator, with respect to a particular ADAPT state.

# Parameters
- `ansatz`: the ADAPT state
- `ADAPT`: the ADAPT protocol. Different protocols may have different scoring strategies.
- `generator`: the generator to be scored
- `H`: the object defining the cost-function
- `ψ0`: an initial quantum state which the `ansatz` operates on

# Returns
- `score`: a scalar number, whose type is `typeof_score(ADAPT)`

# Implementation

In addition to implementing this method (which is mandatory),
    strongly consider over-riding `calculate_scores` also,
    to take advantage of compact measurement protocols,
    or simply the fact that you should only need to evolve your reference state once.

"""
calculate_score(
    ::AbstractAnsatz,
    ::AdaptProtocol,
    ::Generator,
    ::Observable,
    ::QuantumState,
) = NotImplementedError()

"""
    calculate_scores(
        ansatz::AbstractAnsatz,
        ADAPT::AdaptProtocol,
        pool::GeneratorList,
        observable::Observable,
        reference::QuantumState,
    )

Calculate a vector of scores for all generators in the pool.

# Parameters
- `ansatz`: the ADAPT state
- `ADAPT`: the ADAPT protocol. Different protocols may have different scoring strategies.
- `pool`: the list of generators to be scored
- `H`: the object defining the cost-function
- `ψ0`: an initial quantum state which the `ansatz` operates on

# Returns
- `scores`: a vector whose elements are of type `typeof_score(ADAPT)`.
    `scores[i]` is the score for the generator `pool[i]`

"""
function calculate_scores(
    ansatz::AbstractAnsatz,
    ADAPT::AdaptProtocol,
    pool::GeneratorList,
    observable::Observable,
    reference::QuantumState,
)
    scores = Vector{typeof_score(ADAPT)}(undef, length(pool))
    for i in eachindex(pool)
        scores[i] = calculate_score(ansatz, ADAPT, pool[i], observable, reference)
    end
    return scores
end


"""
    adapt!(
        ::AbstractAnsatz,
        ::Trace,
        ::AdaptProtocol,
        ::GeneratorList,
        ::Observable,
        ::QuantumState,
        ::CallbackList,
    )

Update an ansatz with a new generator(s) according to a given ADAPT protocol.

Typically, each call to this function will select a single generator
        whose score has the largest magnitude,
    but richer variants of ADAPT will have richer behavior.

For example, an implementation of Tetris ADAPT would add multiple generators,
    based on both the score and the "disjointness" of the generators
    (which is something that implementation would have to define).

# Parameters
- `ansatz`: the ADAPT state
- `trace`: a history of the ADAPT run thus far
- `ADAPT`: the ADAPT protocol
- `pool`: the list of generators to consider adding to the ansatz
- `H`: the object defining the cost-function
- `ψ0`: an initial quantum state which the `ansatz` operates on
- `callbacks`: a list of functions to be called just prior to updating the ansatz

# Returns
- a `Bool`, indicating whether or not an adaptation was made

# Implementation

Any implementation of this method must be careful to obey the following contract:

1. If your ADAPT protocol decides the ansatz is already converged,
        call `set_converged!(ansatz, true)` and return `false`,
        without calling any callbacks.

2. Fill up a `data` dict with the expensive calculations you have to make anyway.
   See the default implementation for a minimal selection of data to include.

3. BEFORE you actually update the ansatz,
        call each callback in succession, passing it the `data` dict.
   If any callback returns `true`, return `false` without calling any more callbacks,
        and without updating the `ansatz`.

4. After all callbacks have been completed,
    update the ansatz, call `set_optimized!(ansatz, false)`, and return `true`.

Standard operating procedure is to let callbacks do all the updates to `trace`.
Thus, implementations of this method should normally ignore `trace` entirely
    (except in passing it along to the callbacks).
That said, this rule is a "style" guideline, not a contract.

"""
adapt!(
    ::AbstractAnsatz,
    ::Trace,
    ::AdaptProtocol,
    ::GeneratorList,
    ::Observable,
    ::QuantumState,
    ::CallbackList,
) = NotImplementedError()

"""
    (::AbstractCallback)(
        ::Data,
        ::AbstractAnsatz,
        ::Trace,
        ::AdaptProtocol,
        ::GeneratorList,
        ::Observable,
        ::QuantumState,
    )

Callback for adapt iterations, called immediately prior to the ansatz update.

Note that the ansatz is already updated in the optimization callback,
    but not in the adaptation callback.

# Parameters
- Almost all parameters for the `adapt!` method. See that method for details.
- `data`: (replaces `callbacks`) additional calculations the ADAPT method has made
    Keys depend on the protocol. See the `Callbacks` module for some standard choices.

# Returns
- `true` iff ADAPT should terminate, without updating ansatz

"""
(::AbstractCallback)(
    ::Data,
    ::AbstractAnsatz,
    ::Trace,
    ::AdaptProtocol,
    ::GeneratorList,
    ::Observable,
    ::QuantumState,
) = false


##########################################################################################
#= Optimization protocol. =#

"""
    OptimizationProtocol

A distinctive protocol for refining parameters in an ansatz.

# Implementation

Sub-types must implement the following method:
- `optimize!(
    ::AbstractAnsatz,
    ::Trace,
    ::OptimizationProtocol,
    ::Observable,
    ::QuantumState,
    ::CallbackList,
  )::Bool`

"""
abstract type OptimizationProtocol end

"""
    optimize!(
        ansatz::AbstractAnsatz,
        trace::Trace,
        VQE::OptimizationProtocol,
        H::Observable,
        ψ0::QuantumState,
        callbacks::CallbackList,
    )

Update the parameters of an ansatz according to a given optimization protocol.

# Parameters
- `ansatz`: the ADAPT state
- `trace`: a history of the ADAPT run thus far
- `VQE`: the optimization protocol (it doesn't have to be a VQE ^_^)
- `H`: the object defining the cost-function
- `ψ0`: an initial quantum state which the `ansatz` operates on
- `callbacks`: a list of functions to be called just prior to updating the ansatz

# Implementation

Callbacks must be called in each "iteration".
The optimization protocol is free to decide what an "iteration" is,
    but it should generally correspond to "any time the ansatz is changed".
That's not a hard-fast rule, though -
    for example, it doesn't necessarily make sense to call the callbacks
    for each function evaluation in a linesearch.

Any implementation of this method must be careful to obey the following contract:

1. In each iteration, update the ansatz parameters and
        do whatever calculations you need to do.
   Fill up a `data` dict with as much information as possible.
   See the `Callbacks` module for some standard choices.

2. Call each callback in succession, passing it the `data` dict.
   If any callback returns `true`, terminate without calling any more callbacks,
        and discontinue the optimization.

3. After calling all callbacks, check if the ansatz has been flagged as optimized.
   If so, discontinue the optimization.

3. If the optimization protocol terminates successfully without interruption by callbacks,
        call `set_optimized!(ansatz, true)`.
   Be careful to ensure the ansatz parameters actually are the ones found by the optimizer!

Standard operating procedure is to let callbacks do all the updates to `trace`.
Thus, implementations of this method should normally ignore `trace` entirely
    (except in passing it along to the callbacks).
That said, this rule is a "style" guideline, not a contract.

The return type of this method is intentionally unspecified,
    so that implementations can return something helpful for debugging,
    eg. an `Optim` result object.
If the callbacks interrupt your optimization,
    it may be worthwhile to check if they flagged the `ansatz` as converged,
    and modify this return object accordingly if possible.

"""
optimize!(
    ::AbstractAnsatz,
    ::Trace,
    ::OptimizationProtocol,
    ::Observable,
    ::QuantumState,
    ::CallbackList,
) = NotImplementedError()

"""
    (::AbstractCallback)(
        ::Data,
        ::AbstractAnsatz,
        ::Trace,
        ::OptimizationProtocol,
        ::Observable,
        ::QuantumState,
    )

Callback for optimization iterations, called AFTER ansatz update.

Note that the ansatz is already updated in optimization callback,
    but not in the adaptation callback.

# Parameters
- Almost all parameters for the `optimize!` method. See that method for details.
- `data`: (replaces `callbacks`) additional calculations the optimization method has made
    Keys depend on the protocol. See the `Callbacks` module for some standard choices.

# Returns
- `true` iff optimization should terminate

"""
(::AbstractCallback)(
    ::Data,
    ::AbstractAnsatz,
    ::Trace,
    ::OptimizationProtocol,
    ::Observable,
    ::QuantumState,
) = false



##########################################################################################
#= The main entree. =#

"""
    run!(
        ansatz::AbstractAnsatz,
        trace::Trace,
        ADAPT::AdaptProtocol,
        VQE::OptimizationProtocol,
        pool::GeneratorList,
        H::Observable,
        ψ0::QuantumState,
        callbacks::CallbackList,
    )

Loop between optimization and adaptation until convergence.

The `ansatz` and `trace` are mutated throughout,
    so that if the user runs this method in a REPL,
    she can terminate it (eg. by `Ctrl+C`) after however long,
    and still have meaningful results.

# Parameters
- `ansatz`: the ADAPT state
- `trace`: a history of the ADAPT run thus far
- `ADAPT`: the ADAPT protocol
- `VQE`: the optimization protocol (it doesn't have to be a VQE ^_^)
- `pool`: the list of generators to consider adding to the ansatz
- `H`: the object defining the cost-function
- `ψ0`: an initial quantum state which the `ansatz` operates on
- `callbacks`: a list of functions to be called just prior to updating the ansatz

# Returns
- `true` iff the ansatz is converged, with respect to the given protocols and callbacks

"""
function run!(
    ansatz::AbstractAnsatz,
    trace::Trace,
    ADAPT::AdaptProtocol,
    VQE::OptimizationProtocol,
    pool::GeneratorList,
    observable::Observable,
    reference::QuantumState,
    callbacks::CallbackList,
)
    is_converged(ansatz) && return true

    is_optimized(ansatz) || optimize!(ansatz,trace, VQE, observable,reference, callbacks)
    is_optimized(ansatz) || return false

    adapted = adapt!(ansatz, trace, ADAPT, pool, observable, reference, callbacks)
    adapted || return is_converged(ansatz)

    return run!(ansatz, trace, ADAPT, VQE, pool, observable, reference, callbacks)
end