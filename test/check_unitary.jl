#= Simple script to make sure we've done the unitary evolutions right.

Just setup a simple ansatz and a reference, evolve it,
    then convert the ansatz to a unitary and evolve with that,
    and check they give the same final result. =#

import ADAPT

n = 4
N = 1 << n
ψREF = zeros(ComplexF64, N); ψREF[13] = 1

pool = [
    ADAPT.Operators.qubitexcitation(n, 1, 2),
    ADAPT.Operators.qubitexcitation(n, 3, 4),
    ADAPT.Operators.qubitexcitation(n, 1, 2, 3, 4),
]

ansatz = ADAPT.Ansatz(Float64, pool)
push!(ansatz, pool[1] =>  1.5)
push!(ansatz, pool[3] => -0.5)
push!(ansatz, pool[2] =>  0.7)

ψ = ADAPT.evolve_state(ansatz, ψREF)
U = Matrix(N, ansatz)
ψx = U * ψREF

@assert ψ ≈ ψx