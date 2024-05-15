#= Run ADAPT on the Hubbard model with the qubit-excitation pool. =#
#= TODO: Where shall scripts like this live? Not in `test`, I don't think... =#

import ADAPT
import PauliOperators: ScaledPauliVector, FixedPhasePauli, KetBitString, SparseKetBasis

L = 2
u = 0.25

# BUILD OUT THE PROBLEM HAMILTONIAN: a periodic 1d Hubbard lattice at u=0.25
U = 4*u         # Dimensionless parameter u ≡ U/4|t|, and we'll set units so |t|=1.
H = ADAPT.Hamiltonians.hubbard_hamiltonian(L, U, -1.0, pbc=true)

# BUILD OUT THE QUBIT-EXCITATION POOL
pool, target_and_source = ADAPT.Pools.qubitexcitationpool(2L)

# CONSTRUCT A REFERENCE STATE
neel = "0110"^(L >> 1)
(L & 1 == 1) && (neel *= "01")
ket1 = KetBitString{2L}(parse(Int128, neel, base=2))
ψ0 = SparseKetBasis{2L,ComplexF64}(ket1 => 1)

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
    # ADAPT.Callbacks.Serializer(
    #     ansatz_file="dat/ansatz", trace_file="dat/trace",
    #     on_adapt=true, on_iterate=true,
    # ),
    ADAPT.Callbacks.ScoreStopper(1e-3),
    ADAPT.Callbacks.ParameterStopper(10),
]

# RUN THE ALGORITHM
ADAPT.run!(ansatz, trace, adapt, vqe, pool, H, ψ0, callbacks)