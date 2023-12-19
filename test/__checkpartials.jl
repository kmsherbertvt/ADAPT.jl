#= Run ADAPT on the Hubbard model with the qubit-excitation pool. =#
#= TODO: Where shall scripts like this live? Not in `test`, I don't think... =#

import ADAPT
import PauliOperators: ScaledPauliVector, FixedPhasePauli, KetBitString, SparseKetBasis

L = 2
u = 0.25

# BUILD OUT THE PROBLEM HAMILTONIAN: a periodic 1d Hubbard lattice at u=0.25
U = 4*u         # Dimensionless parameter u ≡ U/4|t|, and we'll set units so |t|=1.
H = ADAPT.Operators.hubbard_hamiltonian(L, U, -1.0, pbc=true)

# BUILD OUT THE QUBIT-EXCITATION POOL
pool = ScaledPauliVector{2L}[]
for i in   1:2L
for j in i+1:2L
    push!(pool, ADAPT.Operators.qubitexcitation(2L, i, j))

    for k in j+1:2L
    for l in k+1:2L
        push!(pool, ADAPT.Operators.qubitexcitation(2L, i, j, k, l))
    end; end
end; end

# CONSTRUCT A REFERENCE STATE
#= This reference state is an anti-bonding mixture of the two anti basis states.
    I have chosen this one for this example because it is the first one I found
        for which the qubit excitation pool could start.
=#
neel1 = "0110"^(L >> 1)
(L & 1 == 1) && (neel1 *= "01")
ket1 = KetBitString{2L}(parse(Int128, neel1, base=2))
neel2 = "1001"^(L >> 1)
(L & 1 == 1) && (neel2 *= "10")
ket2 = KetBitString{2L}(parse(Int128, neel2, base=2))
ψ0 = SparseKetBasis{2L,ComplexF64}(ket1 => 1)
# ψ0 = SparseKetBasis{2L,ComplexF64}(ket1 => 1/√2, ket2 => -1/√2)
# ψ0 = zeros(ComplexF64, 1<<2L); ψ0[1+ket1.v] = 1/√2; ψ0[1+ket2.v] = -1/√2

# EVOLVE BY A BIT
evolver = pool[1] + pool[3]
    # A bit arbitrary; simply the generator for which we first noticed an error,
    #                     with ψ0 ~ 0110 - 1001.
display(ψ0)
display(evolver)
ADAPT.evolve_state!(evolver, -1.613, ψ0)
display(ψ0)
println("\n"^4)

# INITIALIZE AN ARBITRARY ANSATZ
ansatz = ADAPT.Ansatz(Float64, pool)
# ansatz = ADAPT.Ansatz(Float64, [pool[3][1]])
# push!(ansatz, pool[3][1] => π/2)
push!(ansatz, pool[1] + pool[3] => -1.613)
push!(ansatz, pool[2] + pool[4] =>  0.237)
push!(ansatz, pool[1] + pool[5] => -0.373)

# MEASURE GRADIENT ANALYTICALLY
@time g0 = ADAPT.gradient(ansatz, H, ψ0)
gO = zero(g0)
fillgO = i -> gO[i] = ADAPT.partial(i, ansatz, H, ψ0)
@time foreach(fillgO, eachindex(ansatz))

# MEASURE GRADIENT NUMERICALLY
import FiniteDifferences
cfd = FiniteDifferences.central_fdm(5, 1)
costfn = ADAPT.Basics.make_costfunction(ansatz, H, ψ0)
gx = FiniteDifferences.grad(cfd, costfn, ADAPT.angles(ansatz))[1]