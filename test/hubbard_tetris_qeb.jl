#= Run TETRIS-ADAPT on the Hubbard model with the qubit-excitation pool. =#

import ADAPT
import PauliOperators: Pauli, PauliSum, ScaledPauli, ScaledPauliVector, FixedPhasePauli, KetBitString, SparseKetBasis

# SYSTEM PARAMETERS
L = 4; u = 0.25

# BUILD OUT THE PROBLEM HAMILTONIAN: an open 1d Hubbard lattice
U = 4*u         # Dimensionless parameter u ≡ U/4|t|, and we'll set units so |t|=1.
H = ADAPT.Hamiltonians.hubbard_hamiltonian(L, U, -1.0, pbc=true)

# EXACT SOLUTION
module Exact
    import ..H, ..L, ..ADAPT
    using LinearAlgebra
    Hm = Matrix(H); E, U_eigs = eigen(Hm) 
    ψ0 = U_eigs[:,1]; E0 = ADAPT.evaluate(H, ψ0);
    println("Exact gs energy = ", E0)
end

# BUILD OUT THE QUBIT-EXCITATION POOL
pool, target_and_source = ADAPT.Pools.qubitexcitationpool(2L)
println("pool size: ", length(pool))

# CONSTRUCT A REFERENCE STATE
neel = "0110"^(L >> 1); (L & 1 == 1) && (neel *= "01")
ket = KetBitString{2L}(parse(Int128, neel, base=2))
ψREF = zeros(ComplexF64, 1<<(2L)); ψREF[1+ket.v] = 1.0

# INITIALIZE THE ANSATZ AND TRACE
ansatz = ADAPT.Ansatz(Float64, pool)
trace = ADAPT.Trace()

# SELECT THE PROTOCOLS
gradient_cutoff = 1e-4; maxiters = 500
adapt = ADAPT.TETRIS_ADAPT.TETRISADAPT(gradient_cutoff)
vqe = ADAPT.OptimOptimizer(:BFGS; g_tol=1e-6)

# SELECT THE CALLBACKS
callbacks = [
    ADAPT.Callbacks.Tracer(:energy, :selected_index, :selected_generator, :selected_score, :scores),
    ADAPT.Callbacks.ParameterTracer(),    
    ADAPT.Callbacks.Printer(:energy, :selected_generator, :selected_score),
    ADAPT.Callbacks.ScoreStopper(gradient_cutoff),
    ADAPT.Callbacks.ParameterStopper(maxiters),
]

# RUN THE ALGORITHM
ADAPT.run!(ansatz, trace, adapt, vqe, pool, H, ψREF, callbacks)

# RESULTS
ψ0 = ADAPT.evolve_state(ansatz, ψREF)
E0 = ADAPT.evaluate(H, ψ0); rel_energy_err = abs((Exact.E0 - E0)/(Exact.E0))
println("relative energy error = ",rel_energy_err)