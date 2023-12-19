#= Finite difference has uncovered an edge case where partial is wrong.
    I presently suspect it is evove_state! that is wrong.
    Test carefully. =#
#= TODO: Where shall scripts like this live? Not in `test`, I don't think... =#

import ADAPT
import PauliOperators: ScaledPauli, FixedPhasePauli, Pauli
import PauliOperators: KetBitString, SparseKetBasis

ψ0 = SparseKetBasis(1, T=ComplexF64); sum!(ψ0, KetBitString(1, 0))
ψ1 = SparseKetBasis(1, T=ComplexF64); sum!(ψ1, KetBitString(1, 1))

# CAREFULLY TEST SINGLE-QUBIT ROTATIONS WHERE WE CAN SEE THE RIGHT ANSWER AT A GLANCE
X = Pauli(1, X=1)
Y = Pauli(1, Y=1)
Z = Pauli(1, Z=1)

display(ADAPT.evolve_state(X, π/2, ψ0))
display(ADAPT.evolve_state(Y, π/2, ψ0))
display(ADAPT.evolve_state(Z, π/2, ψ0))
display(ADAPT.evolve_state(X, π/2, ψ1))
display(ADAPT.evolve_state(Y, π/2, ψ1))
display(ADAPT.evolve_state(Z, π/2, ψ1))
display("")

# SAME BUT FOR SCALED PAULIS
cX = ScaledPauli(X)
cY = ScaledPauli(Y)
cZ = ScaledPauli(Z)

display(ADAPT.evolve_state(cX, π/2, ψ0))
display(ADAPT.evolve_state(cY, π/2, ψ0))
display(ADAPT.evolve_state(cZ, π/2, ψ0))
display(ADAPT.evolve_state(cX, π/2, ψ1))
display(ADAPT.evolve_state(cY, π/2, ψ1))
display(ADAPT.evolve_state(cZ, π/2, ψ1))
display("")

# SAME BUT FOR FIXED-PHASE PAULIS
fX = FixedPhasePauli(1, X=1)
fY = FixedPhasePauli(1, Y=1)
fZ = FixedPhasePauli(1, Z=1)

display(ADAPT.evolve_state(fX, π/2, ψ0))
display(ADAPT.evolve_state(fY, π/2, ψ0))       # NON-HERMITIAN!!!
display(ADAPT.evolve_state(fZ, π/2, ψ0))
display(ADAPT.evolve_state(fX, π/2, ψ1))
display(ADAPT.evolve_state(fY, π/2, ψ1))       # NON-HERMITIAN!!!
display(ADAPT.evolve_state(fZ, π/2, ψ1))
display("")