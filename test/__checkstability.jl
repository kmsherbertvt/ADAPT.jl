#= We've elsewhere uncovered some type instability,
    which may or may not be too blame for a rather egregious-looking memory leak
    whenever we are optimizing something.
I've got to be able to fiddle around with the code here
    to see what it takes to make it go away,
    so let's try to reproduce it here. =#

using Revise
using ADAPT

import PauliOperators: Pauli, ScaledPauli, ScaledPauliVector, PauliSum
import PauliOperators: SparseKetBasis, KetBitString

################################################################################
#= Define work objects. =#

N = 8

adapt = ADAPT.VANILLA
vqe = ADAPT.OptimOptimizer(:BFGS, g_tol=1e-3)
pool = [
    [[1.0 * Pauli(N; Y=q)] for q in 1:N];
    [[1.0 * Pauli(N; X=[q,q+1]), 1.0 * Pauli(N; Y=[q,q+1])] for q in 1:N-1];
]

import Random
Random.seed!(0)
O = let
    H = PauliSum(N)
    for q in 1:N-1
        sum!(H, randn() * Pauli(N; X=q, Y=q+1))
        sum!(H, randn() * Pauli(N; Y=q, Z=q+1))
        sum!(H, randn() * Pauli(N; Z=q, X=q+1))
    end
    H
end

ψREF = let
    ψ = zeros(ComplexF64, 1<<N)
    ψ[13] = 1;
    ψ
end

callbacks = [
    ADAPT.Callbacks.Tracer(
        :energy, :g_norm, :elapsed_time, :elapsed_f_calls, :elapsed_g_calls,
        :selected_index, :selected_score, :scores,
    ),
    ADAPT.Callbacks.ParameterTracer(),
    # ADAPT.Callbacks.Printer(:energy, :g_norm, :selected_index, :selected_score),
    # ADAPT.Callbacks.Serializer(
    #     trace_file="$(vars.outdir)/trace",
    #     on_adapt=true,
    #     on_iterate=true,
    # ),
    ADAPT.Callbacks.ScoreStopper(1e-3),
    ADAPT.Callbacks.ParameterStopper(100),
]

ansatz = ADAPT.Ansatz(Float64, pool)
# push!(ansatz, pool[1] => 1.6)
trace = ADAPT.Trace()

################################################################################
#= Replicate the type instability. =#

@code_warntype ADAPT.evaluate(O, ψREF)
#= It was here.
    Fixed by replacing `PauliOperators.expectation_value`
        with `MyPauliOperators.expectation_value`.
    Shouldn't be a problem at all once we update PauliOperators!
=#




f = ADAPT.make_costfunction(ansatz, O, ψREF)
g! = ADAPT.make_gradfunction!(ansatz, O, ψREF)

x = ADAPT.angles(ansatz)
∇f = similar(x)

@code_warntype f(x)
@code_warntype g!(∇f, x)
#= f(x) had an issue too but fixing evaluate fixed it. =#



state = copy(ψREF)
G = pool[1][1]
θ = 1.2
@code_warntype ADAPT.evolve_state!(G, θ, state)


################################################################################
#= Replicate the apparent memory leak. =#

@time ADAPT.run!(
    ansatz,
    trace,
    adapt,
    vqe,
    pool,
    O,
    ψREF,
    callbacks,
)

#= Results:
With the `evaluate` type instability:   27.97 s, 57.64 M allocations, 40.545 GiB
Without:                                28.52 s, 57.63 M allocations, 40.545 GiB
Using only `MyPauliOperators`:          24.67 s, 13.40 M allocations, 39.007 GiB

Fixing the type instability helps a lot with allocations, but not with total memory.

TODO: I think the problem has _got_ to be from making costates for the "sine branch" and for the gradient evaluation. If you want to solve this, you'll need to make those in-place somehow. Which is not an easy thing. :(

=#