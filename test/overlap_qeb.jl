#= Run Overlap ADAPT on an arbitrary fixed-weight target with the qubit-excitation pool. =#
#= TODO: Where shall scripts like this live? Not in `test`, I don't think... =#

import ADAPT
import PauliOperators: ScaledPauliVector, FixedPhasePauli, KetBitString, SparseKetBasis
import PauliOperators

n = 4       # TOTAL NUMBER OF QUBITS
k = 2       # HAMMING WEIGHT

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

# CONSTRUCT AN ARBITRARY STATEVECTOR FOR TESTING PURPOSES
import Random; Random.seed!(1111)
weight(z) = z == 0 ? 0 : (z & 1) + weight(z>>1)
Φ = SparseKetBasis{n,ComplexF64}()
for z in 0:(1<<n)-1
    (weight(z) != k) && continue                    # ADD ONLY KETS WITH WEIGHT k
    Φ[KetBitString{n}(z)] = randn(ComplexF64)       # RANDOM COMPLEX NUMBER
    # Φ[KetBitString{n}(z)] = randn(Float64)          # RANDOM REAL NUMBER
end
import LinearAlgebra
PauliOperators.scale!(Φ, 1/sqrt(LinearAlgebra.dot(Φ,Φ)))    # NORMALIZE
overlap = ADAPT.OverlapADAPT.Infidelity(Φ)

# CONSTRUCT A REFERENCE STATE
ketstring = "1"^k * "0"^(n-k)
ket = KetBitString{n}(parse(Int128, ketstring, base=2))
ψ0 = SparseKetBasis{n,ComplexF64}(ket => 1)

# INITIALIZE THE ANSATZ AND TRACE
ansatz = ADAPT.Ansatz(Float64, pool)
trace = ADAPT.Trace()

# SELECT THE PROTOCOLS
adapt = ADAPT.VANILLA
vqe = ADAPT.OptimOptimizer(:BFGS; g_tol=1e-6)

# SELECT THE CALLBACKS
callbacks = [
    ADAPT.Callbacks.Tracer(:energy, :selected_index, :selected_score, :scores),
    ADAPT.Callbacks.ParameterTracer(),
    ADAPT.Callbacks.Printer(:energy, :selected_generator, :selected_score),
    ADAPT.Callbacks.ScoreStopper(1e-3),
    ADAPT.Callbacks.ParameterStopper(10),
]

# RUN THE ALGORITHM
ADAPT.run!(ansatz, trace, adapt, vqe, pool, overlap, ψ0, callbacks)
