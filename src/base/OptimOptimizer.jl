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
struct OptimOptimizer
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
    method_kwargs = Dict(Symbol, Any)
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
    cb = options.store_trace ? callback : (trace) -> callback(last(trace))
    kwargs = Dict(field=>getfield(options, field) for field in fieldnames(Optim.Options))
    return Optim.Options(; kwargs..., callback=cb)
end

function make_costfunction(
    ansatz::ADAPT.AbstractAnsatz,
    observable::ADAPT.Observable,
    reference::ADAPT.QuantumState,
    counter::Ref{Int},
)
    return (x) -> (
        counter[] += 1;
        ADAPT.bind!(ansatz, x);
        ADAPT.evaluate(ansatz, observable, reference)
    )
end

function make_gradfunction!(
    ansatz::ADAPT.AbstractAnsatz,
    observable::ADAPT.Observable,
    reference::ADAPT.QuantumState,
    counter::Ref{Int},
)
    return (∇f, x) -> (
        counter[] += 1;
        ADAPT.bind!(ansatz, x);
        for i in eachindex(ansatz);
            ∇f[i] = ADAPT.partial(i, ansatz, observable, reference)
        end;
        ∇f
    )
end

function make_gradfunction(
    ansatz::ADAPT.AbstractAnsatz,
    observable::ADAPT.Observable,
    reference::ADAPT.QuantumState,
    counter::Ref{Int},
)
    g! = make_gradfunction(ansatz, observable, reference, counter)
    return (x) -> (
        ∇f = Vector{ADAPT.energy_type(observable)}(undef, length(x));
        g!(∇f, x);
        ∇f
    )
end

function make_data(state, f_gounter, g_counter)
    return Dict(
        :energy => state.value,
        :g_norm => state.g_norm,
        :elapsed_iterations => state.iteration,
        :elapsed_time => state.metadata["time"],
        :elapsed_f_calls => f_counter[],
        :elapsed_g_calls => g_counter[],
    )
end

function make_callback(
    ansatz::ADAPT.AbstractAnsatz,
    trace::ADAPT.Trace,
    VQE::OptimOptimizer,
    observable::ADAPT.Observable,
    reference::ADAPT.QuantumState,
    callbacks::ADAPT.CallbackList,
    f_counter::Ref{Int},
    g_counter::Ref{Int},
)
    return (state) -> (
        data = make_data(state, f_counter, g_counter);
        stop = false;
        for callback in callbacks;
            stop = stop || callback(data, ansatz, trace, VQE, observable, reference);
            # As soon as `stop` is true, subsequent callbacks are short-circuited.
        end;
        stop || ADAPT.is_optimized(ansatz)
    )

    return store_trace ? (history) -> state_callback(last(history)) : state_callback
end

function ADAPT.optimize!(
    ansatz::ADAPT.AbstractAnsatz,
    trace::ADAPT.Trace,
    VQE::OptimOptimizer,
    observable::ADAPT.Observable,
    reference::ADAPT.QuantumState,
    callbacks::ADAPT.CallbackList,
)
    f_counter = Ref(0)
    g_counter = Ref(0)
    callback = make_callback(
        ansatz, trace, VQE, observable, reference, callbacks,
        f_counter, g_counter,
    )

    result = Optim.optimize(
        make_costfunction(ansatz, observable, reference, f_counter),
        make_gradfunction!(ansatz, observable, reference, g_gounter),
        vqe.method,
        options_with_callback(VQE.options, callback),
    )
    #= TODO: this interface looks different for different optimizers, no?

    You'll have to add a switch on VQE.method <: Optim.FirstOrder..., etc.
    =#

    # IF A CALLBACK DECLARED OPTIMIZATION TO BE SUCCESSFUL, LET THE OPTIMIZER KNOW
    ADAPT.is_optimized(ansatz) && (result.x_converged = true)

    # IF OPTIMIZATION WAS SUCCESSFUL, FLAG `ansatz` AS OPTIMIZED
    Optim.converged(result) && ADAPT.set_optimized!(ansatz, true)

    return result
end