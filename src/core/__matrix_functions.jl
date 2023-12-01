#= DISCLAIMER: This file is meant to be included in ADAPT.jl.
    Names not defined in this file or in base Julia can be found there.
=#

import LinearAlgebra: I

"""
    evolve_state(
        ansatz::AbstractAnsatz,
        unitary::AbstractMatrix{<:Complex},
    )

Calculate the matrix extending the given unitary by an ansatz (on the left).

"""
function evolve_unitary(
    ansatz::AbstractAnsatz,
    unitary::AbstractMatrix{<:Complex},
)
    rotated = deepcopy(unitary)
    evolve_unitary!(ansatz, unitary)
    return rotated
end

"""
    evolve_unitary(
        G::Generator,
        θ::Parameter,
        unitary::AbstractMatrix{<:Complex},
    )

Calculate the matrix extending the given unitary by a single generator (on the left).

"""
function evolve_unitary(
    G::Generator,
    θ::Parameter,
    unitary::AbstractMatrix{<:Complex},
)
    rotated = deepcopy(unitary)
    evolve_unitary!(G, θ, unitary)
    return rotated
end

"""
    evolve_unitary!(
        ansatz::AbstractAnsatz,
        unitary::AbstractMatrix{<:Complex},
    )

Extend a unitary by applying each generator in the ansatz (on the left).

"""
function evolve_unitary!(
    ansatz::AbstractAnsatz,
    unitary::AbstractMatrix{<:Complex},
)
    for (generator, parameter) in ansatz
        evolve_unitary!(generator, parameter, unitary)
    end
    return unitary
end

"""
    evolve_unitary!(
        G::Generator,
        θ::Parameter,
        unitary::AbstractMatrix{<:Complex},
    )

Extend a unitary by applying a single generator (on the left).

"""
function evolve_unitary!(
    G::Generator,
    θ::Parameter,
    unitary::AbstractMatrix{<:Complex},
)
    for i in axes(unitary, 2)
        evolve_state!(G, θ, @view unitary[:,i])
    end
    return unitary
end

"""
    Base.Matrix([F,] N::Int, ansatz::AbstractAnsatz)

Construct the unitary matrix representation of the action of an ansatz.

# Parameters
- `F`: float type; the resulting matrix will be of type Matrix{Complex{F}}
- `N`: size of Hilbert space (ie. the number of rows in the matrix)
- `ansatz`: the ansatz to be represented

"""
function Base.Matrix(F, N::Int, ansatz::AbstractAnsatz)
    unitary = Matrix{Complex{real(eltype(F))}}(I, N, N)
    evolve_unitary!(ansatz, unitary)
    return unitary
end

Base.Matrix(N::Int, ansatz::AbstractAnsatz) = Matrix(Float64, N, ansatz)

"""
    Base.Matrix([F,] N::Int, G::Generator, θ::Parameter)

Construct the unitary matrix representation of ``exp(-iθG)``.

# Parameters
- `F`: float type; the resulting matrix will be of type Matrix{Complex{F}}
- `N`: size of Hilbert space (ie. the number of rows in the matrix)
- `G`: the generator of the matrix
- `θ`: a scalar coefficient multiplying the generator

"""
function Base.Matrix(F, N::Int, G::Generator, θ::Parameter)
    unitary = Matrix{Complex{eltype(F)}}(I, N, N)
    evolve_unitary!(G, θ, unitary)
    return unitary
end

Base.Matrix(N::Int, G::Generator, θ::Parameter) = Matrix(Float64, N, G, θ)