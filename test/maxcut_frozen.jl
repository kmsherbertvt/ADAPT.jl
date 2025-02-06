#= Run ADAPT-QAOA on a MaxCut Hamiltonian,
    so the optimizer only ever has two parameters.

Rather than try to write a new ADAPT ansatz that only optimizes the last two parameters,
    we can do it all with a modified script,
    runing single iterations of ADAPT over and over, with an updated reference state.
The tricky part is stitching the traces together in a seamless way.

=#

import Graphs
import ADAPT

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
prob = 0.5
#
g = Graphs.erdos_renyi(n, prob)
e_list = ADAPT.Hamiltonians.get_unweighted_maxcut(g)

# BUILD OUT THE PROBLEM HAMILTONIAN
H_spv = ADAPT.Hamiltonians.maxcut_hamiltonian(n, e_list)
# Wrap in a QAOAObservable view.
H = ADAPT.ADAPT_QAOA.QAOAObservable(H_spv)

# EXACT DIAGONALIZATION (for accuracy assessment)
#= NOTE: This block now scales as well as evolving the ansatz,
    so there is no need to comment it out. =#
module Exact
    import PauliOperators
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

# SELECT THE CALLBACKS
callbacks = [
    ADAPT.Callbacks.Tracer(:energy, :selected_index, :selected_score, :scores, :g_norm),
    # ADAPT.Callbacks.ParameterTracer(),            # <- Not compatible with this workflow.
    ModalSampleTracer(),
    ADAPT.Callbacks.Printer(:energy, :selected_index, :selected_score),
    ADAPT.Callbacks.ParameterStopper(1),            # <- Essential for this workflow.
]

# BUILD OUT THE STATIC ADAPT OBJECTS
adapt = ADAPT.VANILLA # Can be changed to `ADAPT.Degenerate_ADAPT.DEG_ADAPT`
vqe = ADAPT.OptimOptimizer(:BFGS; g_tol=1e-6, iterations=1000)
pool = ADAPT.ADAPT_QAOA.QAOApools.qaoa_double_pool(n)

# SELECT THE REFERENCE STATE
ψ0 = ones(ComplexF64, 2^n) / sqrt(2^n)

# INITIALIZE THE OBJECTS TO BE UPDATED IN EACH ADAPT ITERATION
γ0 = Ref(0.1)
trace = ADAPT.Trace()
ansatz = ADAPT.ADAPT_QAOA.DiagonalQAOAAnsatz(γ0[], pool, H)
reference = deepcopy(ψ0)

# STOPPING CRITERIA
maxadapt = 100
Gtol = 1e-3
err = 0.5

# RUN THE LOOP
stoppingcriteria = zeros(Bool, 4)
while true
    # INITIALIZE ANSATZ FOR THIS SINGLE ITERATION
    thisansatz = ADAPT.ADAPT_QAOA.DiagonalQAOAAnsatz(γ0[], pool, H)
    if length(ansatz) > 0
        push!(thisansatz.generators, pool[last(trace[:selected_index])])
        push!(thisansatz.β_parameters, zero(γ0[]))  # INITIALIZE β TO ZERO
        push!(thisansatz.γ_parameters, γ0[])        # INITIALIZE γ TO γ0
        ADAPT.set_optimized!(thisansatz, false)
    end

    # RUN ADAPT FOR ONE ITERATION
    success = ADAPT.run!(thisansatz, trace, adapt, vqe, pool, H, reference, callbacks)

    # REGISTER THE RESULTS OF THIS OPTIMIZATION IN OUR RUNNING ANSATZ
    push!(ansatz.generators,   only(thisansatz.generators))
    push!(ansatz.β_parameters, only(thisansatz.β_parameters))
    push!(ansatz.γ_parameters, only(thisansatz.γ_parameters))
    ADAPT.set_optimized!(thisansatz, ADAPT.is_optimized(thisansatz))

    # UPDATE THE OTHER MUTABLES
    reference .= ADAPT.evolve_state(ansatz, ψ0)
    # γ0[] = last(ansatz.γ_parameters)    # OPTIONAL

    # CHECK FOR STOPPING CRITERIA
    stoppingcriteria .= [
        !success,                                   # Optimization failed.
        length(ansatz) ≥ 2maxadapt,                 # Ansatz is too long.
        last(trace[:selected_score]) < Gtol,        # Operators fall below threshold.
        abs(last(trace[:energy]) - Exact.E0) < err  # Close enough to good solution.
        # NOTE: The slow stopper is a bit tedious so I'm not bothering for now...
    ]
    ADAPT.set_converged!(ansatz, any(stoppingcriteria[[3,4]]))
    any(stoppingcriteria) && break
end










# # DISPLAY MOST LIKELY BITSTRINGS FROM EACH ADAPTATION
# println("Most likely bitstrings after each adaption:")
# for z in trace[:modalsample]
#     println(bitstring(z)[end-n+1:end])
# end