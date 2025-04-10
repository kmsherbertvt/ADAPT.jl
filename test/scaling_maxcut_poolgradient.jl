#= Run a single adaptation, and time it, for increasing system sizes. =#

import Graphs
import ADAPT
import LinearAlgebra: norm

import Random; Random.seed!(0)

timings = []

for n in 1:14
    # CREATE THE GRAPH (however you'd like)
    prob = 0.5
    g = Graphs.erdos_renyi(n, prob)
    e_list = ADAPT.Hamiltonians.get_unweighted_maxcut(g)

    # BUILD OUT THE PROBLEM HAMILTONIAN
    H_spv = ADAPT.Hamiltonians.maxcut_hamiltonian(n, e_list)
    # Wrap in a QAOAObservable view.
    H = ADAPT.ADAPT_QAOA.QAOAObservable(H_spv)

    # BUILD OUT THE POOL (however you'd like)
    pool = ADAPT.ADAPT_QAOA.QAOApools.qaoa_double_pool(n)

    # CONSTRUCT A REFERENCE STATE
    ψ0 = ones(ComplexF64, 2^n) / sqrt(2^n); ψ0 /= norm(ψ0)

    # SELECT THE ANSATZ AND ADAPT PROTOCOL
    γ0 = 0.1
    ansatz = ADAPT.ADAPT_QAOA.DiagonalQAOAAnsatz(γ0, pool, H)
    adapt = ADAPT.Degenerate_ADAPT.DegenerateADAPT(1e-8)

    # CALCULATE THE FIRST POOL GRADIENT
    push!(timings, @timed ADAPT.calculate_scores(ansatz, adapt, pool, H, ψ0))

    # YOU CAN ALSO RUN THE WHOLE PROTOCOL FOR THE FIRST ADAPTATION, IF YOU'D PREFER
    #=
    trace = ADAPT.Trace()
    callbacks = [   # Put in whatever you normally do here, to get an accurate timing.
        ADAPT.Callbacks.Tracer(:energy, :selected_index, :selected_score, :scores),
        ADAPT.Callbacks.ParameterTracer(),
        # ModalSampleTracer(),
        ADAPT.Callbacks.Printer(:energy, :selected_index, :selected_score),
        ADAPT.Callbacks.ScoreStopper(1e-3),
        ADAPT.Callbacks.ParameterStopper(100),
        # ADAPT.Callbacks.FloorStopper(0.5, Exact.E0),
        ADAPT.Callbacks.SlowStopper(1.0, 3),
    ]
    push!(timings, @timed ADAPT.adapt!(ansatz, trace, adapt, pool, H, ψ0, callbacks))
    =#

end

# `timings` contains various benchmarking data. See for example:
display([timing.time for timing in timings])