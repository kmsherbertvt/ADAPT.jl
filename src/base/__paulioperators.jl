################################################################################
#= TODO: These belong in PauliOperators, probably. =#
################################################################################

using PauliOperators
using PauliOperators: AbstractPauli
AnyPauli = Union{AbstractPauli, PauliSum, ScaledPauliVector}

import LinearAlgebra

"""
    Sigh... this seems to not be defined for any AbstractPauli.
    Best do so. For now all I need is FixedPhasePauli, to get the durned thing to run.
    There's more, though. FixedPhasePauli * Vector ?? I'm confused. Giving up for now...
"""
function LinearAlgebra.mul!(C::Vector{T}, A::FixedPhasePauli{N}, B::Vector{T}) where {T,N}
    ndim = size(B,1)

    # Check dimensions
    size(B) == size(C) || throw(DimensionMismatch)
    ndim == 2^N || throw(DimensionMismatch)

    op = A; coeff = 1
    #= @inbounds @simd =# for i in 0:ndim-1
        (phase, j) = op * KetBitString{N}(i)
        C[j.v+1] += phase * coeff * B[i+1]
    end
end

""" Of course this one is missing... ^_^
Note strict typing in out, because Paulis themselves are strictly typed.
"""
function Base.:*(ps::ScaledPauliVector{N}, ψ::SparseKetBasis{N,T}) where {N,T}
    σ = SparseKetBasis(N, T=ComplexF64)
    for p in ps
        sum!(σ, p*ψ)
    end
    return σ
end


function Base.sum!(ψ1::SparseKetBasis{N,T}, ψ2::SparseKetBasis{N,T}) where {N,T}
    for (ket, c) in ψ2
        sum!(ψ1, ket, c)
    end
end

function PauliOperators.clip!(ψ::SparseKetBasis{N,T}; thresh=1e-16) where {N,T}
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

function PauliOperators.expectation_value(
    pauli::AnyPauli,
    state::Union{SparseKetBasis,AbstractVector},
)
    bra = pauli * state
    #= TODO: Allocations! Should use `mul!(tmp, pauli, state)` somehow... =#
    return LinearAlgebra.dot(state, bra)
end

function braket(
    pauli::AnyPauli,
    bra::Union{SparseKetBasis,AbstractVector},
    ket::Union{SparseKetBasis,AbstractVector},
)
    covector = pauli * ket
    #= TODO: Allocations! Should use `mul!(tmp, pauli, ket)` somehow... =#
    return LinearAlgebra.dot(bra, covector)
end

#= TODO: lmul!. Probably a separate method for each Pauli type and each quantum state. =#


################################################################################
#= TODO: These are heavier-weight and may or may not be appropriate. =#
################################################################################

"""
Cross-type multiplication.
Best to discourage ever doing this operation.
Needed for a lazy commutator, but not necessarily needed long-term.
We'll return pauli sum for now.
"""
function Base.:*(ps1::PauliSum{N}, ps2::ScaledPauliVector{N}) where {N}
    out = PauliSum(N)
    for (op1, coeff1) in ps1.ops

        for scaledpauli in ps2
            op2 = scaledpauli.pauli
            coeff2 = scaledpauli.coeff

            prod = op1 * op2
            if haskey(out, prod)
                out[prod] += get_phase(op1, op2)*coeff1*coeff2
            else
                out[prod] = get_phase(op1, op2)*coeff1*coeff2
            end
            # out.ops[prod] = get(out.ops, prod) + get_phase(prod)*coeff1*coeff2
        end
    end
    return out
end
function Base.:*(ps1::ScaledPauliVector{N}, ps2::PauliSum{N}) where {N}
    out = PauliSum(N)
    for scaledpauli in ps1
        op1 = scaledpauli.pauli
        coeff1 = scaledpauli.coeff

        for (op2, coeff2) in ps2.ops
            prod = op1 * op2
            if haskey(out, prod)
                out[prod] += get_phase(op1, op2)*coeff1*coeff2
            else
                out[prod] = get_phase(op1, op2)*coeff1*coeff2
            end
            # out.ops[prod] = get(out.ops, prod) + get_phase(prod)*coeff1*coeff2
        end
    end
    return out
end

function Base.:*(ps1::PauliSum{N}, op2::FixedPhasePauli{N}) where {N}
    out = PauliSum(N)
    for (op1, coeff1) in ps1.ops
        prod = op1 * op2
        if haskey(out, prod)
            out[prod] += get_phase(op1, op2)*coeff1
        else
            out[prod] = get_phase(op1, op2)*coeff1
        end
    end
    return out
end
function Base.:*(op1::FixedPhasePauli{N}, ps2::PauliSum{N}) where {N}
    out = PauliSum(N)
    for (op2, coeff2) in ps2.ops
        prod = op1 * op2
        if haskey(out, prod)
            out[prod] += get_phase(op1, op2)*coeff2
        else
            out[prod] = get_phase(op1, op2)*coeff2
        end
    end
    return out
end

"""
    measure_commutator(
        A::AnyPauli,
        B::AnyPauli,
        Ψ::Union{SparseKetBasis,AbstractVector},
    )

Calculate the expectation value of the commutator, ie. ⟨Ψ|[A,B]|Ψ⟩.

TODO: There *could* be a place for this in PauliOperators,
        but it would need to be carefully fleshed out type by type.
    A and B needn't be Hermitian in general (though I assume they are here),
        so my intuition is rather lacking.

"""
function measure_commutator(
    A::AnyPauli,
    B::AnyPauli,
    Ψ::Union{SparseKetBasis,AbstractVector},
)
    commutator = A*B - B*A
    #= TODO: ALLOCATIONS!
        Diksha has a much more efficient version in the ACSE package.
        But, it would need distinct methods for all the type combinations.
        All my operators are fixed size,
            so I don't know that I care about the allocations for ADAPT.jl.
        But I think the more efficient implementations pry do belong in PauliOperators.jl.
    =#
    return expectation_value(commutator, Ψ)
end







#= TODO: Calculate DLA

1. basis(pool): make a vector of all the Paulis appearing in pool
2. as_matrix(pool, basis): make a matrix of each operator in pool, in representation of basis (most useful as input to `rank`)
3. generate_dynamicalliealgebra(pool): loop through using below
4. fill_dynamicalliealgebra!(algebra, queue, pool, basis): deal with one element of pool, updating algebra, queue, and pool. Then return all args for convenient looping

=#