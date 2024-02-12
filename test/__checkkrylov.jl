#= Working out syntax for using KrylovKit.exponentiate with PauliSum. =#

using PauliOperators
using LinearAlgebra
using KrylovKit

# Yo ho, yo ho.
# Base.getindex(k::SparseKetBasis{N,T}, i::Integer) where {N,T} = get(k.coeffs, KetBitString{N}(i-1), zero(T))
Base.similar(k::SparseKetBasis) = zero(k)
LinearAlgebra.rmul!(k::SparseKetBasis, c::Number) = (PauliOperators.scale!(k, c); k)
    # TODO: scale! ought to return k already. (But KrylovKit *needs* rmul! too.)
# TODO: Final version should constrain N type parameter to be the same for these:
LinearAlgebra.mul!(w::SparseKetBasis, v::SparseKetBasis, α::Number) = (
    empty!(w.coeffs);
    for (ket, c) in v.coeffs;
        w[ket] = α*c;
    end;
    w
)
LinearAlgebra.axpby!(α, v::SparseKetBasis, β, w::SparseKetBasis) = (
    for (ket, c) in v.coeffs;
        w[ket] = β*get(w, ket, 0) + α * c
    end;
    w
)
LinearAlgebra.axpy!(α, v::SparseKetBasis, w::SparseKetBasis) = LinearAlgebra.axpby!(α, v, 1, w)

##########
#=
- Make a PauliSum G of non-commuting ScaledPaulis.
- Make a random low-weight SparseKetBasis x.

- Convert G into a matrix.
- Convert x into a vector.
- Evalute G⋅x and exp(-itG)⋅x by brute force linear algebra.

- Wrapper function to evaluate G⋅x, for arbitrary statevector x.
- Call function for x, to match against G⋅x linear algebra.
- Call KrylovKit.exponentiate(A, -it, x), to match against exp(-itG)⋅x linear algebra.

- Wrapper function to evaluate G⋅x, for arbitrary ket x.
- Call function for x then convert to vector, to match against G⋅x linear algebra.
- Call KrylovKit.exponentiate(A, -it, x) then convert to vector, to match against ...
=#

n = 6       # INCREASE THIS TO SEE SCALING, I GUESS

# CONSTRUCT PAULI SUM ``aXX + bYY`` ON THREE QUBITS
G = PauliSum(n)
sum!(G,  .34 * Pauli(n; X=[1, 2]))
sum!(G, -.89 * Pauli(n; X=1, Y=2))

# CONSTRUCT INITIAL STATEVECTOR ``a|100⟩ + b|010⟩ + c|001⟩``
ψ = SparseKetBasis(n; T=ComplexF64);
ψ[KetBitString(n, 4)] = 0.62 + .15im
ψ[KetBitString(n, 2)] = 0.18 - .32im
ψ[KetBitString(n, 1)] = 0.55 + .00im
PauliOperators.scale!(ψ, 1/norm(values(ψ)))

# SELECT ROTATION ANGLE
t = 0.76

# "KNOWN SOLUTION" LINEAR ALGEBRA
G_ = Matrix(G)
ψ_ = Vector(ψ)
Gψ_ = G_ * ψ_
expGψ_ = cis(-t .* G_) * ψ_

# LINEAR FUNCTION
A(x) = G * x

# TRY IT OUT FOR VECTORS
@assert A(ψ_) ≈ Gψ_
@assert exponentiate(A, -im*t, ψ_)[1] ≈ expGψ_

# TRY IT OUT FOR SPARSE KET BASES
@assert Vector(A(ψ)) ≈ Gψ_
@assert Vector(exponentiate(A, -im*t, ψ)[1]) ≈ expGψ_
# NOTE: Too bad, SparseKetBasis does not sufficiently "behave as a vector". Not worth it.

# TEST THE PACKAGE'S EVOLUTION METHOD
import ADAPT
state = copy(ψ_)
@assert ADAPT.evolve_state(G, t, state) ≈ expGψ_
ADAPT.evolve_state!(G, t, state)
@assert state ≈ expGψ_