#= DISCLAIMER: This file is meant to be included in ADAPT.jl.
    Names not defined in this file or in base Julia can be found there.
=#

"""
    Trace

Semantic alias for a compact record of the entire ADAPT run.

Keys can in principle be any Symbol at all.
You can design your own protocols to fill the data,
    and your own callbacks to use it.
That said, see the `Callbacks` module for some standard choices.

"""
Trace = Dict{Symbol, Any}

"""
    Data

Semantic alias for trace-worthy information from a single adapt or vqe iteration.

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

- Printers: print some data to the standard output stream
- Tracers: save some data to the running trace
- Stoppers: check some data to see if we've hit a stopping condition

In particular, the standard way to converge an adapt run is to include a `ScoreStopper`.
Otherwise, the run will keep adding parameters until every score is essentially zero.

More details can be found in the `Callbacks` module,
    where many standard callbacks are already implemented.

# Implementation

Callbacks are implemented as callable objects, with two choices of method header
    (one for adaptations, one for optimization iterations).

- `(::AbstractCallback)(
    ::AbstractAnsatz,
    ::Trace,
    ::AdaptProtocol,
    ::GeneratorList,
    ::Observable,
    ::QuantumState,
    ::Data,
  )`

- `(::AbstractCallback)(
    ::AbstractAnsatz,
    ::Trace,
    ::OptimizationProtocol,
    ::Observable,
    ::QuantumState,
    ::Data,
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

"""
    (::AbstractCallback)(
        ::AbstractAnsatz,
        ::Trace,
        ::AdaptProtocol,
        ::GeneratorList,
        ::Observable,
        ::QuantumState,
        ::Data,
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
    ::AbstractAnsatz,
    ::Trace,
    ::AdaptProtocol,
    ::GeneratorList,
    ::Observable,
    ::QuantumState,
    ::Data,
) = false

"""
    (::AbstractCallback)(
        ::AbstractAnsatz,
        ::Trace,
        ::OptimizationProtocol,
        ::Observable,
        ::QuantumState,
        ::Data,
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
    ::AbstractAnsatz,
    ::Trace,
    ::OptimizationProtocol,
    ::Observable,
    ::QuantumState,
    ::Data,
) = false



