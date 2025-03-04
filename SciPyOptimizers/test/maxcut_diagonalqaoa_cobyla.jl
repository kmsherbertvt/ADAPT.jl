#= Run ADAPT-QAOA on a MaxCut Hamiltonian.

This script is identical to `maxcut_diagonalqaoa.jl`
    except that it uses a `SciPyOptimizer` to run COBYLA,
    instead of the `OptimOptimizer` running BFGS.

=#

import Graphs
import ADAPT
import SciPyOptimizers
import PauliOperators
import PauliOperators: ScaledPauliVector, FixedPhasePauli, KetBitString, SparseKetBasis
import PauliOperators: PauliSum
import LinearAlgebra: norm

import Random; Random.seed!(0)

"""
    ModalSampleTracer()

At each adaptation, identify the most likely bitstring and save it as an integer.

In the context of QAOA, this identifies the most reasonable partition.

"""
struct ModalSampleTracer <: ADAPT.AbstractCallback end

function (tracer::ModalSampleTracer)(
    ::ADAPT.Data, ansatz::ADAPT.AbstractAnsatz, trace::ADAPT.Trace,
    ::ADAPT.AdaptProtocol, ::ADAPT.GeneratorList,
    ::ADAPT.Observable, ψ0::ADAPT.QuantumState,
)
    ψ = ADAPT.evolve_state(ansatz, ψ0)      # THE FINAL STATEVECTOR
    imode = argmax(abs2.(ψ))                # MOST LIKELY INDEX
    zmode = imode-1                         # MOST LIKELY BITSTRING (as int)

    push!( get!(trace, :modalsample, Any[]), zmode )
    return false
end

# DEFINE A GRAPH
n = 12

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
#= NOTE: This block now scales as well as evolving the ansatz,
    so there is no need to comment it out. =#
module Exact
    import ..PauliOperators
    import ..H, ..n
    Emin = Ref(Inf); ketmin = Ref(PauliOperators.KetBitString{n}(0))
    for v in 0:1<<n-1
        ket = PauliOperators.KetBitString{n}(v)
        vec = PauliOperators.SparseKetBasis{n,ComplexF64}(ket => 1)
        Ev = real((H*vec)[ket])
        if Ev < Emin[]
            Emin[] = Ev
            ketmin[] = ket
        end
    end
    ψ0 = Vector(PauliOperators.SparseKetBasis{n,ComplexF64}(ketmin[] => 1))
    E0 = Emin[]

    ρ = abs2.(ψ0)                           # THE FINAL PROBABILITY DISTRIBUTION
    pmax, imax = findmax(ρ)
    ketmax = PauliOperators.KetBitString(n, imax-1) # THE MOST LIKELY BITSTRING
end
println("Exact ground-state energy: ",Exact.E0)
println("Best cut: ",Exact.ketmax)

# BUILD OUT THE POOL
pool = ADAPT.ADAPT_QAOA.QAOApools.qaoa_double_pool(n)

# ANOTHER POOL OPTION
#pool = ADAPT.Pools.two_local_pool(n)

println("Generator data type: ", typeof(pool[1]))
println("Note: using the 'diagonal' ADAPT-QAOA ansatz, the observable type is fixed.")

# CONSTRUCT A REFERENCE STATE
ψ0 = ones(ComplexF64, 2^n) / sqrt(2^n); ψ0 /= norm(ψ0)

# INITIALIZE THE ANSATZ AND TRACE
ansatz = ADAPT.ADAPT_QAOA.DiagonalQAOAAnsatz(0.1, pool, H)
# the first argument (a hyperparameter) can in principle be set to values other than 0.1
trace = ADAPT.Trace()

# SELECT THE PROTOCOLS
adapt = ADAPT.Degenerate_ADAPT.DegenerateADAPT(1e-8)  # Can be changed to `ADAPT.VANILLA`
vqe = SciPyOptimizers.SciPyOptimizer("COBYLA";
    tol=1e-6,       # Not exactly sure what this entails, honestly...
    rhobeg=2π/10,   # Reasonable initial changes to values.
                    #= NOTE: I made this one (2π/10) up.

                    Maybe it's a good choice for proper angles like β; I dunno.
                    It's not necessarily appropriate for γ and β
                        to use the same value here,
                        but I can't help that without implementing COBYLA myself. =#
    maxiter=1000,   # Maximum number of iterations per optimization.
    disp=true,      # Use scipy's own verbose printing.
)

# SELECT THE CALLBACKS
callbacks = [
    ADAPT.Callbacks.Tracer(:energy, :selected_index, :selected_score, :scores),
    ADAPT.Callbacks.ParameterTracer(),
    ModalSampleTracer(),
    ADAPT.Callbacks.Printer(:energy, :selected_index, :selected_score),
    ADAPT.Callbacks.ScoreStopper(1e-3),
    ADAPT.Callbacks.ParameterStopper(100),
    ADAPT.Callbacks.FloorStopper(0.5, Exact.E0),
    ADAPT.Callbacks.SlowStopper(1.0, 3),
]

# RUN THE ALGORITHM
success = ADAPT.run!(ansatz, trace, adapt, vqe, pool, H, ψ0, callbacks)
println(success ? "Success!" : "Failure - optimization didn't converge.")

# DISPLAY MOST LIKELY BITSTRINGS FROM EACH ADAPTATION
println("Most likely bitstrings after each adaption:")
for z in trace[:modalsample]
    println(bitstring(z)[end-n+1:end])
end