import ..ADAPT

import Optim

"""
    OptimOptimizer(method, options)

# Parameters
- `method`: an optimizer object from the `Optim` package
- `options`: an options object from the `Optim` package

IMPORTANT: the `callback` attribute of `options` will be generated automatically
        whenever `ADAPT.optimize!` is called, to insert dynamic callbacks.
    If you provide your own callback in `options`, it will be ignored.
    Use the `ADAPT.Callback` framework to gain extra behavior throughout optimization.
    If this framework does not meet your needs,
        you'll need to implement your own `OptimizationProtocol`.

"""
struct OptimOptimizer <: ADAPT.OptimizationProtocol
    method::Optim.AbstractOptimizer
    options::Optim.Options
end

"""
    OptimOptimizer(method::Symbol; options...)

A convenience constructor to create `OptimOptimizers` without referring to `Optim`.

# Parameters
- `method`: a symbol-ization of the `Optim` method

# Kewyord Arguments
You can pass any keyword argument accepted either by
    your `Optim` method's constructor, or by that of `Optim.Options`.
If you try to pass a `callback` keyword argument, it will be ignored (see above).

"""
function OptimOptimizer(method::Symbol; options...)
    method_type = getfield(Optim, method)

    # EXTRACT METHOD kwargs FROM options
    method_fields = fieldnames(method_type)
    method_kwargs = Dict{Symbol, Any}()
    for field in keys(options)
        if field in method_fields
            method_fields[field] = pop!(options, field)
        end
    end

    # CONSTRUCT THE `Optim` OBJECTS
    return OptimOptimizer(
        method_type(; method_kwargs...),
        Optim.Options(; options...),
    )
end

function options_with_callback(options, callback)
    cb = options.store_trace ? (trace) -> callback(last(trace)) : callback
    kwargs = Dict(field=>getfield(options, field) for field in fieldnames(Optim.Options))
    return Optim.Options(; kwargs..., callback=cb)
end

function make_costfunction(
    ansatz::ADAPT.AbstractAnsatz,
    observable::ADAPT.Observable,
    reference::ADAPT.QuantumState,
)
    x0 = copy(ADAPT.angles(ansatz))
    return (x) -> (
        x0 .= ADAPT.angles(ansatz);         # SAVE THE ORIGINAL STATE
        ADAPT.bind!(ansatz, x);
        f = ADAPT.evaluate(ansatz, observable, reference);
        ADAPT.bind!(ansatz, x0);            # RESTORE THE ORIGINAL STATE
        f
    )
end

function make_gradfunction!(
    ansatz::ADAPT.AbstractAnsatz,
    observable::ADAPT.Observable,
    reference::ADAPT.QuantumState,
)
    x0 = copy(ADAPT.angles(ansatz))
    return (∇f, x) -> (
        x0 .= ADAPT.angles(ansatz);         # SAVE THE ORIGINAL STATE
        ADAPT.bind!(ansatz, x);
        ADAPT.gradient!(∇f, ansatz, observable, reference);
        ADAPT.bind!(ansatz, x0);            # RESTORE THE ORIGINAL STATE
        ∇f
    )
end

function make_gradfunction(
    ansatz::ADAPT.AbstractAnsatz,
    observable::ADAPT.Observable,
    reference::ADAPT.QuantumState,
    counter::Ref{Int},
)
    g! = make_gradfunction!(ansatz, observable, reference)
    return (x) -> (
        ∇f = Vector{ADAPT.typeof_energy(observable)}(undef, length(x));
        g!(∇f, x);
        ∇f
    )
end

function make_data(iterdata, objective, state)
    return ADAPT.Data(
        :energy => iterdata.value,
        :g_norm => iterdata.g_norm,
        :elapsed_iterations => iterdata.iteration,
        :elapsed_time => iterdata.metadata["time"],
        :elapsed_f_calls => only(objective.f_calls),
        :elapsed_g_calls => only(objective.df_calls),
    )
end

function make_callback(
    ansatz::ADAPT.AbstractAnsatz,
    trace::ADAPT.Trace,
    VQE::OptimOptimizer,
    observable::ADAPT.Observable,
    reference::ADAPT.QuantumState,
    callbacks::ADAPT.CallbackList,
    objective::Optim.OnceDifferentiable,
    state::Optim.AbstractOptimizerState,
)
    iterdata_callback = (iterdata) -> (
        ADAPT.bind!(ansatz, state.x);   # TODO: Is `x` a guaranteed field of `state`?
        data = make_data(iterdata, objective, state);
        stop = false;
        for callback in callbacks;
            stop = stop || callback(data, ansatz, trace, VQE, observable, reference);
            # As soon as `stop` is true, subsequent callbacks are short-circuited.
        end;
        # TODO: Consider updating `state.x` so callbacks can control x. Seems dangerous...
        stop || ADAPT.is_optimized(ansatz)
    )

    history_callback = (history) -> iterdata_callback(last(history))

    return VQE.options.store_trace ? history_callback : iterdata_callback
end

function ADAPT.optimize!(
    ansatz::ADAPT.AbstractAnsatz,
    trace::ADAPT.Trace,
    VQE::OptimOptimizer,
    observable::ADAPT.Observable,
    reference::ADAPT.QuantumState,
    callbacks::ADAPT.CallbackList,
)
    # INITIALIZE OPTIMIZATION OBJECTS
    objective = Optim.OnceDifferentiable(
        make_costfunction(ansatz, observable, reference),
        make_gradfunction!(ansatz, observable, reference),
        copy(ADAPT.angles(ansatz)),     # NOTE: This copy is pry redundant. But safer.
    )
    # TODO: Generalize interface for 0th/2nd order methods, and 1st w/finite difference.

    state = Optim.initial_state(
        VQE.method,
        VQE.options,                    # NOTE: Pretty sure this argument is inert.
        objective,
        copy(ADAPT.angles(ansatz)),
    )

    callback = make_callback(
        ansatz, trace, VQE, observable, reference, callbacks,
        objective, state,
    )

    # RUN OPTIMIZATION
    result = Optim.optimize(
        objective,
        copy(ADAPT.angles(ansatz)),
        VQE.method,
        options_with_callback(VQE.options, callback),
        state,
    )

    # IF A CALLBACK DECLARED OPTIMIZATION TO BE SUCCESSFUL, LET THE OPTIMIZER KNOW
    ADAPT.is_optimized(ansatz) && (result.x_converged = true)

    # IF OPTIMIZATION WAS SUCCESSFUL, FLAG `ansatz` AS OPTIMIZED
    Optim.converged(result) && ADAPT.set_optimized!(ansatz, true)

    return result
end