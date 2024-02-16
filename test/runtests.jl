using ADAPT
using Test

import PauliOperators: Pauli, ScaledPauli, ScaledPauliVector, PauliSum
import PauliOperators: SparseKetBasis, KetBitString

##########################################################################################
#= DEFINE RE-USABLE ADAPT OBJECTS =#

N = 4

BFGS = ADAPT.OptimOptimizer(:BFGS, g_tol=1e-3)

pools = (
    PauliSum = [
        [PauliSum(Pauli(N; Y=q)) for q in 1:N];
        [Pauli(N; X=[q,q+1]) + Pauli(N; Y=[q,q+1]) for q in 1:N-1];
    ],
    ScaledPauliVector = [
        [[1.0 * Pauli(N; Y=q)] for q in 1:N];
        [[1.0 * Pauli(N; X=[q,q+1]), 1.0 * Pauli(N; Y=[q,q+1])] for q in 1:N-1];
    ],
    # TODO: ScaledPauli
    # TODO: Pauli
)

observables = (
    PauliSum = let
        H = PauliSum(N)
        for q in 1:N-1
            sum!(H, Pauli(N; X=q))
            sum!(H, Pauli(N; Z=[q,q+1]))
        end
        sum!(H, Pauli(N; X=N))
        H
    end,
    # TODO: ScaledPauliVector
    # TODO: ScaledPauli
    # TODO: Pauli
    Infidelity = let
        ψ = zeros(ComplexF64, 1<<N)
        ψ[4] = 1/√2
        ψ[13] = 1/√2
        ADAPT.OverlapADAPT.Infidelity(ψ)
    end,
)

references = (
    SparseKetBasis = let
        ket = KetBitString{N}(parse(Int128, "1100", base=2))
        sparseket = SparseKetBasis{N,ComplexF64}(ket => 1)
        sparseket
    end,
    Vector = let
        ψ = zeros(ComplexF64, 1<<N)
        ψ[13] = 1;
        ψ
    end,
)

##########################################################################################
#= FUNCTION TO RUN A VALIDATION, FOR A SPECIFIC COMBO =#

function run_tests(combo)
    label = join(map(string, combo), ".")

    ansatztype = ansatze[combo[1]]
    adapt = adapts[combo[2]]
    vqe = vqes[combo[3]]
    pool = pools[combo[4]]
    observable = observables[combo[5]]
    reference = references[combo[6]]

    ansatz = ansatztype(Float64, pool)

    ADAPT.validate(ansatz, adapt, vqe, pool, observable, reference; label=label)
end


##########################################################################################
#= VALIDATE SELECT COMBINATIONS =#

@testset "ADAPT.jl" begin
    @testset "Basics" begin
        ADAPT.validate(
            ADAPT.Ansatz(Float64, pools[:ScaledPauliVector]),
            ADAPT.VANILLA,
            BFGS,
            pools[:ScaledPauliVector],
            observables[:PauliSum],
            references[:Vector];
            label = "ScaledPauli[] Pool, Statevector",
        )

        ADAPT.validate(
            ADAPT.Ansatz(Float64, pools[:ScaledPauliVector]),
            ADAPT.VANILLA,
            BFGS,
            pools[:ScaledPauliVector],
            observables[:PauliSum],
            references[:SparseKetBasis];
            label = "ScaledPauli[] Pool, SparseKetBasis",
        )

        ADAPT.validate(
            ADAPT.Ansatz(Float64, pools[:PauliSum]),
            ADAPT.VANILLA,
            BFGS,
            pools[:PauliSum],
            observables[:PauliSum],
            references[:Vector];
            label = "PauliSum Pool, Statevector",
        )

        # ADAPT.validate(
        #     ADAPT.Ansatz(Float64, pools[:PauliSum]),
        #     ADAPT.VANILLA,
        #     BFGS,
        #     pools[:PauliSum],
        #     observables[:PauliSum],
        #     references[:SparseKetBasis];
        #     label = "PauliSum Pool, SparseKetBasis",
        # )
        # TODO: Just some linear algebra, but SparseKetBasis is changing drastically soon.
    end

    @testset "Optimization-Free" begin
        ADAPT.validate(
            ADAPT.Ansatz(Float64, pools[:ScaledPauliVector]),
            ADAPT.VANILLA,
            ADAPT.OptimizationFreeADAPT.OPTIMIZATION_FREE,
            pools[:ScaledPauliVector],
            observables[:PauliSum],
            references[:Vector];
            label = "ScaledPauli[] Pool, Statevector",
        )
    end

    @testset "Overlap" begin
        overlap = let
            ψ = zeros(ComplexF64, 1<<N)
            ψ[4] = 1/√2
            ψ[13] = 1/√2
            ADAPT.OverlapADAPT.Infidelity(ψ)
        end

        ADAPT.validate(
            ADAPT.Ansatz(Float64, pools[:ScaledPauliVector]),
            ADAPT.VANILLA,
            BFGS,
            pools[:ScaledPauliVector],
            overlap,
            references[:Vector];
            label = "Statevector",
        )

        ADAPT.validate(
            ADAPT.Ansatz(Float64, pools[:ScaledPauliVector]),
            ADAPT.VANILLA,
            BFGS,
            pools[:ScaledPauliVector],
            overlap,
            references[:SparseKetBasis];
            label = "SparseKetBasis",
        )
    end
end
