module SciPyOptimizers
    import ADAPT
    import SciPy

    """
        SciPyOptimizer(method, options)

    An optimization protocol attempting to wrap around `scipy.optimize.minimize`.

    # Parameters
    - `method`: a string, identifying the optimization method (e.g. `BFGS` or `COBYLA`)
    - `bounds`: the bounds to pass for the call to `scipy.optimize.minimize`
    - `constraints`: the constraints to pass for the call to `scipy.optimize.minimize`
    - `tol`: the tolerance to pass for the call to `scipy.optimize.minimize`
    - `options`: the options to pass for the call to `scipy.optimize.minimize`

    """
    struct SciPyOptimizer <: ADAPT.OptimizationProtocol
        method::String
        bounds::Any
        constraints::Any
        tol::Any
        options::Any
    end

    """
        SciPyOptimizer(method; bounds, constraints, tol, options...)

    A convenience constructor to create `SciPyOptimizers` with default arguments.

    The `bounds`, `constraints`, and `tol` kwargs will default to `nothing`,
        so the call to `scipy.optimize.minimize` is equivalent to not passing them at all,
        if you choose to omit them.
    I can offer no advice on how to construct them in a way compatible with `SciPy.jl`
        (especially `constraints`...).

    The remaining kwargs will be wrapped up neatly
        into a `Dict` for the `options` parameter.

    """
    function SciPyOptimizer(
        method::String;
        bounds=nothing,
        constraints=nothing,
        tol=nothing,
        options...,
    )
        return SciPyOptimizer(method, bounds, constraints, tol, Dict(options...))
    end

    #=

    Scipy has two callback signatures,
        a modern one taking an object with `x` and `fun` attributes,
        and a legacy one just takes `x` itself.
    The following functions semantically check which kind of argument it is
        and extracts the information appropriately,
        calling the function with `x` if necessary.

    =#
    current_parameters(arg, fun) = hasproperty(arg, :x) ? arg.x : arg
    current_energy(arg, fun) = hasproperty(arg, :fun) ? arg.fun : fun(arg)

    function make_callback(
        ansatz::ADAPT.AbstractAnsatz,
        trace::ADAPT.Trace,
        VQE::SciPyOptimizer,
        observable::ADAPT.Observable,
        reference::ADAPT.QuantumState,
        callbacks::ADAPT.CallbackList,
        fun::Function,
    )
        return (intermediate_result) -> (
            x = current_parameters(intermediate_result, fun);
            E = current_energy(intermediate_result, fun);
            data = ADAPT.Data(:energy=>E);
            ADAPT.bind!(ansatz, x);
            stop = false;
            for callback in callbacks;
                stop = stop || callback(data, ansatz, trace, VQE, observable, reference);
                # As soon as `stop` is true, subsequent callbacks are short-circuited.
            end;
            stop || ADAPT.is_optimized(ansatz)
            #= TODO: I don't know how to implement callback termination with `SciPy`.
                With the more modern signature,
                    I'm supposed to be able to raise `StopIteration`.
                But I don't know how to do this with `SciPy.jl` or even `PyCall.jl`
                    (but it should be doable somehow!).

                I think it isn't possible at all using the legacy callback,
                    which is what COBYLA always uses. =#
            #= TODO: The modern signature relies on introspection; does it work at all?
                That would mean `SciPy.jl` is feeding everything into `scipy`
                    with unadulterated names; I'm not sure that is possible... =#
        )
    end

    function ADAPT.optimize!(
        ansatz::ADAPT.AbstractAnsatz,
        trace::ADAPT.Trace,
        VQE::SciPyOptimizer,
        observable::ADAPT.Observable,
        reference::ADAPT.QuantumState,
        callbacks::ADAPT.CallbackList,
    )
        # ASSEMBLE OPTIMIZATION OBJECTS
        fun = ADAPT.make_costfunction(ansatz, observable, reference)
        jac = ADAPT.make_gradfunction(ansatz, observable, reference)
        x0 = collect(ADAPT.angles(ansatz))
        callback = make_callback(ansatz, trace, VQE, observable, reference, callbacks, fun)

        # RUN THE OPTIMIZATION
        result = SciPy.optimize.minimize(
            fun, x0,
            # No need for `args`, with this framework.
            method = VQE.method,
            jac = jac,  # NOTE: Not used for gradient-free methods like COBYLA.
            # I don't think we'll ever bother to support `hess` or `hessp`.
            bounds = VQE.bounds,
            constraints = VQE.constraints,
            tol = VQE.tol,
            options = VQE.options,
            callback = callback,
        )

        display(result)

        # IF A CALLBACK DECLARED OPTIMIZATION TO BE SUCCESSFUL, LET THE OPTIMIZER KNOW
        ADAPT.is_optimized(ansatz) && (result["success"] = true)

        # IF OPTIMIZATION WAS SUCCESSFUL, FLAG `ansatz` AS OPTIMIZED
        result["success"] && ADAPT.set_optimized!(ansatz, true)

        return result
    end

end # module SciPyOptimizers
