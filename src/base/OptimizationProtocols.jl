"""
    OptimizationFree

The optimization protocol which just doesn't do anything.

There are no iterations, so there is no reason to callback. Contract obliged!

"""
struct OptimizationFree <: OptimizationProtocol end
OPTIMIZATION_FREE = OptimizationFree()

optimize!(
    ansatz::AbstractAnsatz,
    ::Trace,
    ::OptimizationFree,
    ::Observable,
    ::QuantumState,
    ::CallbackList,
) = set_optimized!(ansatz, true)