import ..ADAPT

"""
    OptimizationFree

The optimization protocol which just doesn't do anything.

There are no iterations, so there is no reason to callback. Contract obliged!

"""
struct OptimizationFree <: ADAPT.OptimizationProtocol end
OPTIMIZATION_FREE = OptimizationFree()

function ADAPT.optimize!(
    ansatz::ADAPT.AbstractAnsatz,
    trace::ADAPT.Trace,
    VQE::OptimizationFree,
    observable::ADAPT.Observable,
    reference::ADAPT.QuantumState,
    callbacks::ADAPT.CallbackList,
)
    # MEASURE ENERGY
    energy = ADAPT.evaluate(ansatz, observable, reference)
    data = ADAPT.Data(:energy => energy)

    # CALL ONE ROUND OF CALLBACKS
    stop = false
    for callback in callbacks
        stop = stop || callback(data, ansatz, trace, VQE, observable, reference)
        # As soon as `stop` is true, subsequent callbacks are short-circuited.
    end

    # TERMINATE OPTIMIZATION
    ADAPT.set_optimized!(ansatz, true)
end