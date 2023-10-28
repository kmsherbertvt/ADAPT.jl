import ..ADAPT

import LinearAlgebra

import PauliOperators
import PauliOperators: KetBitString, SparseKetBasis
import PauliOperators: AbstractPauli
import PauliOperators: FixedPhasePauli, ScaledPauliVector
import PauliOperators: PauliSum, ScaledPauli

################################################################################
#= TODO: These belong in PauliOperators, probably. =#
################################################################################
function Base.sum!(ψ1::SparseKetBasis{N,T}, ψ2::SparseKetBasis{N,T}) where {N,T}
    for (ket, c) in ψ2
        sum!(ψ1, ket, c)
    end
end

function clip!(ψ::SparseKetBasis{N,T}; thresh=1e-16) where {N,T}
    filter!(p -> abs(p.second) > thresh, ψ)
end

function rotate!(
    state::AbstractVector,
    pauli::FixedPhasePauli,
    angle::Number,
)
    sinebranch = pauli * state
    #= TODO: Allocations! Should use `mul!(tmp, G, ψ0)` somehow... =#
    sinebranch .*= sin(angle)
    state .*= cos(angle)
    state .+= sinebranch
end

function rotate!(
    state::SparseKetBasis,
    pauli::FixedPhasePauli,
    angle::Number,
)
    sinebranch = pauli * state
    PauliOperators.scale!(sinebranch, sin(angle))
    PauliOperators.scale!(state, cos(angle))
    sum!(state, sinebranch)
end

function expectation_value(
    pauli::AbstractPauli,
    state::Union{SparseKetBasis,AbstractVector},
)
    bra = pauli * state
    #= TODO: Allocations! Should use `mul!(tmp, pauli, state)` somehow... =#
    return LinearAlgebra.dot(state, bra)
end

function braket(
    pauli::AbstractPauli,
    bra::Union{SparseKetBasis,AbstractVector},
    ket::Union{SparseKetBasis,AbstractVector},
)
    covector = pauli * ket
    #= TODO: Allocations! Should use `mul!(tmp, pauli, ket)` somehow... =#
    return LinearAlgebra.dot(bra, covector)
end

function LinearAlgebra.lmul!(
    pauli::AbstractPauli,
    state::SparseKetBasis,
)
    #= TODO: Alas, I think this needs to be done separately for each Pauli type. =#
end

################################################################################

ADAPT.typeof_energy(::AbstractPauli) = Float64
#= NOTE:
    This type assertion assumes two things:
    1. Scaled AbstractPaulis are hard-coded to used `ComplexF64`.
    2. The observable is Hermitian (so its energy is guaranteed to be real).

    Strictly speaking, using the plain `Pauli` type undermines both assumptions,
        since half of them are anti-Hermitian and they don't have *any* float type.
    So...try to avoid using `Pauli` as an `Observable`?
    It should work fine if the phase is real, but...it'd be less than robust...
=#

function ADAPT.evolve_state!(
    G::FixedPhasePauli,
    θ::ADAPT.Parameter,
    Ψ::ADAPT.QuantumState,
)
    rotate!(Ψ, G, θ)
    Ψ isa SparseKetBasis && clip!(Ψ)
end

function ADAPT.evolve_state!(
    G::ScaledPauliVector,
    θ::ADAPT.Parameter,
    Ψ::ADAPT.QuantumState,
)
    for scaled in G
        rotate!(Ψ, scaled.pauli, scaled.coeff * θ)
    end
end

function ADAPT.evolve_state!(
    G::PauliSum,
    θ::ADAPT.Parameter,
    Ψ::ADAPT.QuantumState,
)
    #= TODO: Implement this, using KrylovKit.exponentiate!

    https://jutho.github.io/KrylovKit.jl/stable/man/matfun/#KrylovKit.exponentiate

    This requires `x` to be "any Julia type with vector like behavior".
    Not too sure *which* behavior is required.
    So, there's a chance `PauliSum` will only be compatible with `AbstractVector`,
        but that's okay.

    =#
    return NotImplementedError()
end



function ADAPT.evaluate(
    H::AbstractPauli,
    Ψ::ADAPT.QuantumState,
)
    return real(expectation_value(H, Ψ))
end


"""
    partial(
        index::Int,
        ansatz::AbstractAnsatz,
        observable::Observable,
        reference::QuantumState,
    )

The partial derivative of a cost-function with respect to the i-th parameter in an ansatz.

The ansatz is assumed to apply a unitary rotation `exp(-iθG)`,
    where `G` is the (Hermitian) generator,
    and generators with a lower index are applied to the state earlier.
Ansatz sub-types may change both behaviors.

# Parameters
- `index`: the index of the parameter to calculate within `ansatz`
- `ansatz`: the ADAPT state
- `H`: the object defining the cost-function
- `ψ0`: an initial quantum state which the `ansatz` operates on

# Returns
- a number of type `typeof_energy(observable)`.

"""
function ADAPT.partial(
    index::Int,
    ansatz::ADAPT.AbstractAnsatz,
    observable::ADAPT.Observable,
    reference::ADAPT.QuantumState,
)
    state = deepcopy(reference)

    # PARTIAL EVOLUTION
    for i in 1:index
        generator, parameter = ansatz[i]
        evolve_state!(generator, parameter, state)
    end

    # REFLECTION
    costate = pauli * state

    # FINISH EVOLUTION
    for i in 1+index:length(ansatz)
        generator, parameter = ansatz[i]
        evolve_state!(generator, parameter, state)
        evolve_state!(generator, parameter, costate)
    end

    return -2 * real(braket(observable, costate, state))
end