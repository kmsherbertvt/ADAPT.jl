#= Check the infidelity gradients against finite difference,
    and check that adding a new operator gives a gradient matching scores. =#

import ADAPT
import PauliOperators: ScaledPauliVector, FixedPhasePauli, KetBitString, SparseKetBasis
import PauliOperators: clip!

n = 4

# BUILD OUT THE QUBIT-EXCITATION POOL
pool = ScaledPauliVector{n}[]
for i in   1:n
for j in i+1:n
    push!(pool, ADAPT.Operators.qubitexcitation(n, i, j))

    for k in j+1:n
    for l in k+1:n
        push!(pool, ADAPT.Operators.qubitexcitation(n, i, j, k, l))
    end; end
end; end
#= NOTE: This pool conserves Hamming weight, so it obv. won't find a random target.
    But it is sufficient for testing purposes.
    A full-fledged example script would
        restrict the random target to have support only on a particular Hamming sector.
    =#

function as_sparseketbasis(ψ::AbstractVector)
    N = length(ψ)
    n = round(Int, log2(N))
    ψ_ = SparseKetBasis{n, eltype(ψ)}()
    for i in eachindex(ψ)
        ψ_[KetBitString{n}(i-1)] = ψ[i]
    end
    clip!(ψ_)
    return ψ_
end


# CONSTRUCT AN ARBITRARY STATEVECTOR FOR TESTING PURPOSES
import Random; Random.seed!(1111)
Φ = randn(ComplexF64, 1<<n)
Φ = as_sparseketbasis(Φ)                        # CAST TO SPARSE STATEVECTOR
overlap = ADAPT.OverlapADAPT.Infidelity(Φ)

# CONSTRUCT A REFERENCE STATE
ψ0 = zeros(ComplexF64, 1<<n); ψ0[2] = 1
ψ0 = as_sparseketbasis(ψ0)                      # CAST TO SPARSE STATEVECTOR
    #= NOTE: This reference belongs to the single-Hamming weight sector,
        since 2 maps to ket 0..01. =#

# INITIALIZE AN ARBITRARY ANSATZ
ansatz = ADAPT.Ansatz(Float64, pool)
push!(ansatz, pool[1] => -1.613)
push!(ansatz, pool[2] =>  0.237)
push!(ansatz, pool[1] => -0.373)

# MEASURE GRADIENT ANALYTICALLY
@time g0 = ADAPT.gradient(ansatz, overlap, ψ0)
gO = zero(g0)
fillgO = i -> gO[i] = ADAPT.partial(i, ansatz, overlap, ψ0)
@time foreach(fillgO, eachindex(ansatz))

# MEASURE GRADIENT NUMERICALLY
import FiniteDifferences
cfd = FiniteDifferences.central_fdm(5, 1)
costfn = ADAPT.Basics.make_costfunction(ansatz, overlap, ψ0)
gx = FiniteDifferences.grad(cfd, costfn, ADAPT.angles(ansatz))[1]

# MEASURE SCORES AS INTENDED
@time s0 = ADAPT.calculate_scores(ansatz, ADAPT.VANILLA, pool, overlap, ψ0)
sO = zero(s0)
fillsO = i -> sO[i] = ADAPT.calculate_score(ansatz, ADAPT.VANILLA, pool[i], overlap, ψ0)
@time foreach(fillsO, eachindex(pool))

# MEASURE SCORES VIA THE GRADIENT
sx = zero(s0)
for i in eachindex(pool)
    candidate = deepcopy(ansatz); push!(candidate, pool[i] => 0.0)
    sx[i] = last(ADAPT.gradient(candidate, overlap, ψ0))
end
