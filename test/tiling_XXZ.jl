#= Run ADAPT on the XXZ model with the Pauli pool. =#

import ADAPT
import PauliOperators: Pauli, PauliSum, ScaledPauli, ScaledPauliVector, ≈

# SYSTEM PARAMETERS
L_small = 3; L_large = 4; Jxy = 1.0; Jz = 0.5; PBCs = false
pooltype="fullpauli"  # fullpauli || qubitadapt || qubitexcitation

# BUILD OUT THE PROBLEM HAMILTONIAN: an open XXZ model
XXZHam = ADAPT.Hamiltonians.xyz_model(L_small, Jxy, Jxy, Jz, PBCs)

# CONSTRUCT A REFERENCE STATE
neel = "01"^(L_small >> 1); (L_small & 1 == 1) && (neel *= "0"); neel_index = parse(Int128, neel, base=2)

# BUILD OUT THE POOL
if pooltype=="fullpauli"
    pool = ADAPT.Pools.fullpauli(L_small)
elseif pooltype == "qubitexcitation"
    pool, target_and_source = ADAPT.Pools.qubitexcitationpool(L_small)
elseif pooltype == "qubitadapt"
    pool = ADAPT.Pools.qubitadaptpool(L_small)
end
println("pool size ",length(pool)); # println("pool: ",pool)

# SELECT THE PROTOCOLS
adapt = ADAPT.Degenerate_ADAPT.DEG_ADAPT
vqe = ADAPT.OptimOptimizer(:BFGS; g_tol=1e-6)
gradient_cutoff = 1e-4; maxiters = 500

# SELECT THE CALLBACKS
callbacks = [
    ADAPT.Callbacks.Tracer(:energy, :selected_generator, :selected_index, :selected_score, :scores),
#     ADAPT.Callbacks.Printer(:energy, :selected_generator, :selected_score),
    ADAPT.Callbacks.ScoreStopper(gradient_cutoff),
    ADAPT.Callbacks.ParameterStopper(maxiters),
]

# EXACT DIAGONALIZATION
module Exact
    import ..XXZHam, ..L_small
    using LinearAlgebra
    Hm = Matrix(XXZHam); E, U = eigen(Hm) # NOTE: Comment out after first run when debugging.
    ψ0 = U[:,1]
    E0 = E[1]
end
println("Exact gs energy: ",Exact.E0)

# RUN MANY ADAPT-VQE TRIALS ON THE SMALL PROBLEM INSTANCE
chosen_operators = ScaledPauliVector{L_small}[]
trials = 100
println("running trials on the small problem instance...")  
for i=1:trials
    # INITIALIZE THE REFERENCE STATE
    ψ0 = zeros(ComplexF64,2^L_small); 
    ψ0[neel_index+1] = 1.0    # Neel

    # INITIALIZE THE ANSATZ AND TRACE
    ansatz = ADAPT.Ansatz(Float64, pool)
    trace = ADAPT.Trace()

    # RUN THE ALGORITHM
    ADAPT.run!(ansatz, trace, adapt, vqe, pool, XXZHam, ψ0, callbacks)
    selected_operators = trace[:selected_generator][1:end-1]
#     energy_err = abs((Exact.E0-(trace[:energy][end-1]))/Exact.E0); print("  ",energy_err)
    for selected_op in selected_operators
        duplicate = false
        for chosen_op in chosen_operators
            if selected_op ≈ chosen_op
                duplicate = true
                break
            end
        end
        if !duplicate
            push!(chosen_operators, selected_op)
        end
    end
end
println("\nTrials on small problem instance complete.")
println(length(chosen_operators), " chosen operators: ",chosen_operators,"\n")

# # FIND DLA OF OPERATORS
# t_0_Lie = time(); Lie_alg = Lie_algebra_elements(chosen_operators); t_f_Lie = time(); dt = (t_f_Lie - t_0_Lie)/60.0
# # println("walltime to calculate Lie algebra for L = $L_tile was $dt minutes\n")
# println("DLA of chosen operators ",length(Lie_alg), " ",Lie_alg,"\n")


# RUN ADAPT-VQE ON THE LARGE PROBLEM INSTANCE

# BUILD OUT THE PROBLEM HAMILTONIAN: an open XXZ model
XXZHam = ADAPT.Hamiltonians.xyz_model(L_large, Jxy, Jxy, Jz, PBCs)

# EXACT DIAGONALIZATION
module Exact_large
    import ..XXZHam, ..L_large
    using LinearAlgebra
    Hm = Matrix(XXZHam); E, U = eigen(Hm) # NOTE: Comment out after first run when debugging.
    ψ0 = U[:,1]
    E0 = E[1]
end
println("Exact gs energy: ",Exact_large.E0)

# CONSTRUCT A REFERENCE STATE
neel = "01"^(L_large >> 1); (L_large & 1 == 1) && (neel *= "0"); neel_index = parse(Int128, neel, base=2)
ψ0 = zeros(ComplexF64,2^L_large); ψ0[neel_index+1] = 1.0

# BUILD OUT THE POOL
# Lie_algebra_elements(chosen_op_strings) # calculate the Lie algebra of the chosen operators
pool = ADAPT.Pools.tile_operators(L_small, L_large, chosen_operators, PBCs)

# SELECT THE PROTOCOLS
adapt = ADAPT.VANILLA
vqe = ADAPT.OptimOptimizer(:BFGS; g_tol=1e-6)

# INITIALIZE THE ANSATZ AND TRACE
ansatz = ADAPT.Ansatz(Float64, pool)
trace = ADAPT.Trace()

# RUN THE ALGORITHM
ADAPT.run!(ansatz, trace, adapt, vqe, pool, XXZHam, ψ0, callbacks)

# RESULTS
ψEND = ADAPT.evolve_state(ansatz, ψ0)
num_ADAPT_iters = length(ansatz)
E0 = trace[:energy][end-1]; rel_energy_err = abs((Exact_large.E0 - E0)/(Exact_large.E0))
println("final ADAPT energy = ", E0)
println("ansatz length: ", length(ansatz))
println("relative energy error = ",rel_energy_err)