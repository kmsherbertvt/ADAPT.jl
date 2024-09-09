#= Run ADAPT-QAOA on a MaxCut Hamiltonian.

This script is identical to `maxcut_qaoa.jl` except that it uses the `DiagonalQAOAAnsatz`,
    which is a more sophisticated (i.e. objectively better) implementation.
At some point when this package is strictly versioned,
    the new implementation will replace the old
    (I'd have done so already but the constructor is not backwards compatible).

=#

import Graphs
import ADAPT
import PauliOperators: ScaledPauliVector, FixedPhasePauli, KetBitString, SparseKetBasis
import PauliOperators: PauliSum
import LinearAlgebra: norm

import Random; Random.seed!(0)

# DEFINE A GRAPH
n = 6

# EXAMPLE OF ERDOS-RENYI GRAPH
prob = 0.5
g = Graphs.erdos_renyi(n, prob)

# EXAMPLE OF ANOTHER ERDOS-RENYI
#ne = 7
#g = Graphs.erdos_renyi(n, ne)

# EXTRACT MAXCUT FROM GRAPH
e_list = ADAPT.Hamiltonians.get_unweighted_maxcut(g)

# BUILD OUT THE PROBLEM HAMILTONIAN
H_spv = ADAPT.Hamiltonians.maxcut_hamiltonian(n, e_list)
# Wrap in a QAOAObservable view.
H = ADAPT.ADAPT_QAOA.QAOAObservable(H_spv)


# ANOTHER WAY TO BUILD OUT THE PROBLEM HAMILTONIAN
#d = 3 # degree of regular graph
#H = ADAPT.Hamiltonians.MaxCut.random_regular_max_cut_hamiltonian(n, d)
println("Observable data type: ",typeof(H))

# EXACT DIAGONALIZATION (for accuracy assessment)
#= NOTE: Comment this out to avoid building exponentially-sized matrices!
    (but then you'll have to edit or remove the FloorStopper callback. =#
module Exact
    import ..H
    using LinearAlgebra
    Hm = Matrix(H); E, U = eigen(Hm) # NOTE: Comment out after first run when debugging.
    ψ0 = U[:,1]
    E0 = real(E[1])
end
println("Exact ground-state energy: ",Exact.E0)

# BUILD OUT THE POOL
pool = ADAPT.ADAPT_QAOA.QAOApools.qaoa_double_pool(n)

# ANOTHER POOL OPTION
#pool = ADAPT.Pools.two_local_pool(n)

println("Generator data type: ", typeof(pool[1]))
println("Note: using the 'diagonal' ADAPT-QAOA ansatz, the observable type is fixed.")

# CONSTRUCT A REFERENCE STATE
ψ0 = ones(ComplexF64, 2^n) / sqrt(2^n); ψ0 /= norm(ψ0)

# INITIALIZE THE ANSATZ AND TRACE
ansatz = ADAPT.ADAPT_QAOA.MixedQAOAAnsatz(0.1, pool, H)
# the first argument (a hyperparameter) can in principle be set to values other than 0.1
trace = ADAPT.Trace()

# SELECT THE PROTOCOLS
adapt = ADAPT.VANILLA # Can be changed to `ADAPT.Degenerate_ADAPT.DEG_ADAPT`
vqe = ADAPT.OptimOptimizer(:BFGS; g_tol=1e-6)
    #= NOTE: Add `iterations=10` to set max iterations per optimization loop. =#

# SELECT THE CALLBACKS
callbacks = [
    ADAPT.Callbacks.Tracer(:energy, :selected_index, :selected_score, :scores),
    ADAPT.Callbacks.ParameterTracer(),
    ADAPT.Callbacks.Printer(:energy, :selected_index, :selected_score),
    ADAPT.Callbacks.ScoreStopper(1e-3),
    ADAPT.Callbacks.ParameterStopper(100),
    ADAPT.Callbacks.FloorStopper(0.5, Exact.E0),
    ADAPT.Callbacks.SlowStopper(1.0, 3),
]

# RUN THE ALGORITHM
success = ADAPT.run!(ansatz, trace, adapt, vqe, pool, H, ψ0, callbacks)
println(success ? "Success!" : "Failure - optimization didn't converge.")
