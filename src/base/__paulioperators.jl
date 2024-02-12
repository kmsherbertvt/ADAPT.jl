################################################################################
#= TODO: These belong in PauliOperators, probably. =#
################################################################################

using PauliOperators
using PauliOperators: AbstractPauli
AnyPauli = Union{AbstractPauli, PauliSum, ScaledPauliVector}

import LinearAlgebra

#=
TODO: Mutating functions (ie. those ending in !) should return the mutated object.
This is a matter of style and conveniente to allow chaining,
    eg. `clip!(sum!(paulisum, pauli))`
    but I've also just discovered it is necessary for correct functionality
    of things like `reduce(sum!, paulis; init=paulisum)`.
=#

#= TODO: Any function args with Matrix or Vector should be AbstractMatrix/Vector.
    Case in point: they should operate perfectly well on views. >_>
=#
Base.:*(p::AnyPauli, in::AbstractArray) = Base.:*(p, Array(in))

function LinearAlgebra.mul!(
    C::AbstractVector{T},
    A::AbstractPauli{N},
    B::AbstractVector{T},
) where {T,N}
    ndim = size(B,1)

    # Check dimensions
    size(B) == size(C) || throw(DimensionMismatch)
    ndim == 2^N || throw(DimensionMismatch)

    C .= 0
    op = A; coeff = 1
    #= @inbounds @simd =# for i in 0:ndim-1
        (phase, j) = op * KetBitString{N}(i)
        C[j.v+1] += phase * coeff * B[i+1]
    end
end

function Base.zero(::SparseKetBasis{N,T}) where {N,T}
    return SparseKetBasis(N, T=T)
end

"""
TODO: Consult with Nick before adding this definition to PauliOperators.

I hesitate for two reasons:

1. It is not "lazy". It allocates a new array.
    Not unprecedented but not ideal.
    Not sure the proper way to make it lazy.

2. Column vector adjoint should properly be a row vector, rather than reversed.
    Can't think of why we'd ever use ScaledPauliVector as a column vector,
        but its data type is so, properly.

But, this definition achieves desired polymorphism in evolving by ScaledPauliVector,
    so if Nick okays it, I'm happy with it.
The alternative is a dedicated `unevolve` function with a tedious special case
    for unevolving ansatze whose generators are ScaledPauliVector...

"""
function Base.adjoint(ps::ScaledPauliVector)
    return [p' for p in reverse(ps)]
end

"""
TODO: This adjoint is not strictly "lazy". But I don't think anyone will care.
"""
function Base.adjoint(ψ::SparseKetBasis)
        ψ_ = zero(ψ)
        for (ket, coeff) in ψ
            ψ_[ket] = coeff'
        end
        return ψ_
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
    for (ket, c) in ψ2.coeffs
        sum!(ψ1, ket, c)
    end
end

function PauliOperators.clip!(ψ::SparseKetBasis{N,T}; thresh=1e-16) where {N,T}
    filter!(p -> abs(p.second) > thresh, ψ.coeffs)
end

"""
TODO: VERY SPECIFICALLY ASSERT that pauli xz=00 is to be interpreted as I,
                                    pauli xz=10 is to be interpreted as X,
                                    pauli xz=01 is to be interpreted as Z,
                                and pauli xz=11 is to be interpreted as Y,
                                despite the last usually being interpreted as iY.
    Also clear this definition with Nick before putting it in his package...
"""
function cis!(state::AbstractVector, pauli::FixedPhasePauli, angle)
    sinebranch = pauli * state
    sinebranch .*= PauliOperators.get_phase(pauli) * im * sin(angle)
    state .*= cos(angle)
    state .+= sinebranch
end

function cis!(state::SparseKetBasis, pauli::FixedPhasePauli, angle)
    sinebranch = pauli * state
    PauliOperators.scale!(sinebranch, PauliOperators.get_phase(pauli) * im * sin(angle))
    PauliOperators.scale!(state, cos(angle))
    sum!(state, sinebranch)
end

function cis!(state::Union{SparseKetBasis,AbstractVector}, pauli::Pauli, angle)
    return cis!(state, pauli.pauli, angle * PauliOperators.get_phase(pauli))
end

function cis!(state::Union{SparseKetBasis,AbstractVector}, pauli::ScaledPauli, angle)
    return cis!(state, pauli.pauli, angle * pauli.coeff)
end

function PauliOperators.expectation_value(
    pauli::AnyPauli,
    state::Union{SparseKetBasis,AbstractVector},
)
    covector = pauli * state
    #= TODO: Allocations! Should use `mul!(tmp, pauli, state)` somehow... =#
    return LinearAlgebra.dot(state, covector)
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