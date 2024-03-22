#= DISCLAIMER: This file is meant to be included in ADAPT.jl.
    Names not defined in this file or in base Julia can be found there.
=#

import LinearAlgebra
import FiniteDifferences

import Test
import Test: @testset, @test

"""
    validate(
        ansatz::AbstractAnsatz,
        adapt::AdaptProtocol,
        vqe::OptimizationProtocol,
        pool::GeneratorList,
        observable::Observable,
        reference::QuantumState;
        kwargs...
    )

Validate that ADAPT will work correctly with the given types.

It actually *runs* ADAPT,
    so ensure your pool and observable are as simple as the types allow.

The mandatory arguments are exactly those found in the `run!` method,
    except there is no trace.

Keyword Arguments
-----------------
- label: the name of the test-set (useful when validating more than one set of types).
- tolerance: the default tolerance for numerical tests
- evolution: special tolerance for the evolution test, or `nothing` to skip
- evaluation: special tolerance for the evaluation test, or `nothing` to skip
- gradient: special tolerance for the gradient test, or `nothing` to skip
- scores: special tolerance for the scores test, or `nothing` to skip

"""
function validate(
    ansatz::AbstractAnsatz,
    adapt::AdaptProtocol,
    vqe::OptimizationProtocol,
    pool::GeneratorList,
    observable::Observable,
    reference::QuantumState;
    label= "ADAPT Validation",
    tolerance=1e-10,
    evolution=tolerance,
    evaluation=tolerance,
    gradient=tolerance,
    scores=tolerance,
)
    @testset "$label" begin

    # VALIDATE RUNTIME
    runtime_tests = validate_runtime(
        ansatz, adapt, vqe, pool, observable, reference;
        verbose=false,
    )

    # CHECK FOR ERRORS IN THE RUNTIME VALIDATION
    # runtime_validated = !runtime_tests.anynonpass
    #= NOTE: By all appearances, the above line should be sufficient.
            But for some reason, the `anynonpass` field is not set correctly.
            I guess this field is set ex-post facto, for some reason?

        The following, manual solution is specific
            to how I've structured the runtime tests.
        If that structure changes, so must this code...
    =#
    runtime_validated = Ref(true)
    for testset in runtime_tests.results
        isa(testset, Test.Error) && (runtime_validated[] = false)
        for resultwrapper in testset.results
            isa(resultwrapper, Test.Error) && (runtime_validated[] = false)
            for result in resultwrapper.results
                isa(result, Test.Error) && (runtime_validated[] = false)
                isa(result, Test.Fail) && (runtime_validated[] = false)
            end
        end
    end

    # IF SUCCESSFUL, CHECK FOR ERRORS OR TEST FAILURES
    if runtime_validated[]
        validate_consistency(ansatz, adapt, pool, observable, reference)

        isnothing(evolution) || validate_evolution(
                pool[1], one(typeof_parameter(ansatz)), reference; tolerance=evolution)
        isnothing(evaluation) || validate_evaluation(
                observable, reference; tolerance=evaluation)
        isnothing(gradient) || validate_gradient(
                ansatz, observable, reference; tolerance=gradient)
        isnothing(scores) || validate_scores(
                ansatz, adapt, pool, observable, reference; tolerance=scores)
    end

    end # @testset
end


"""
    @runtime do_time, ex
    @runtime(do_time, ex)

A macro to check that an expression evaluates without error,
    optionally including an explicit test for runtime.
"""
macro runtime(do_time, ex)
    quote
        @testset $(string(ex)) begin
            @test try $(esc(do_time)) ? @showtime($(esc(ex))) : $(esc(ex))
                true
            catch e
                if isa(e, ADAPT.NotImplementedError)
                    println(stderr, "*** $(e.msg) ***")
                end
                false
            end
        end
    end
end

"""
    validate_runtime(
        ansatz::AbstractAnsatz,
        adapt::AdaptProtocol,
        vqe::OptimizationProtocol,
        pool::GeneratorList,
        observable::Observable,
        reference::QuantumState;
        verbose=true,
    )

Check that every core ADAPT function can run for the given types.

If `verbose` is true, this method also explicilty @time's everything,
    to catch any super-obvious memory leaks when called manually.

Note that this *will* run ADAPT for one iteration,
    so ensure your pool and observable are as simple as the types allow.

"""
function validate_runtime(
    ansatz::AbstractAnsatz,
    adapt::AdaptProtocol,
    vqe::OptimizationProtocol,
    pool::GeneratorList,
    observable::Observable,
    reference::QuantumState;
    verbose=true,
)
    v = verbose     # A shorter alias
    return @testset "Runtime Validation" begin

    # INITIAL SETUP, to ensure we have something to test
    isempty(pool) && error("Cannot validate ADAPT with an empty pool.")
    generator = first(pool)
    angle = 1.0
    push!(ansatz, generator => angle)           # Guarantees evolution *does* something.
    callbacks = [ADAPT.Callbacks.ParameterStopper(1)]   # Immediately terminates ADAPT.

    # RUNTIME TESTS
    v && println("--- Validating Runtime ---")
    v && println("Runtime from `__number_aliases`:")
    @testset "Number Aliases" begin; end

    v && println("Runtime from `__quantum_objects`:")
    @testset "Quantum Objects" begin
        @runtime v typeof_energy(observable)
    end

    v && println("Runtime from `__quantum_functions`:")
    @testset "Quantum Functions" begin
        @runtime v evolve_state(ansatz, reference)
        @runtime v evolve_state(generator, angle, reference)
        @runtime v evolve_state!(ansatz, deepcopy(reference))
        @runtime v evolve_state!(generator, angle, deepcopy(reference))

        @runtime v evaluate(ansatz, observable, reference)
        @runtime v evaluate(observable, reference)

        @runtime v partial(1, ansatz, observable, reference)
        @runtime v gradient(ansatz, observable, reference)
        result = zeros(typeof_energy(observable), length(ansatz))
        @runtime v gradient!(result, ansatz, observable, reference)

        @runtime v make_costfunction(ansatz, observable, reference)
        @runtime v make_gradfunction(ansatz, observable, reference)
        @runtime v make_gradfunction!(ansatz, observable, reference)
    end

    v && println("Runtime from `__ansatz`:")
    @testset "Ansatz" begin
        @runtime v ADAPT.__get__generators(ansatz)
        @runtime v ADAPT.__get__parameters(ansatz)
        @runtime v ADAPT.__get__optimized(ansatz)
        @runtime v ADAPT.__get__converged(ansatz)

        @runtime v typeof_parameter(ansatz)

        @runtime v is_optimized(ansatz)
        @runtime v is_converged(ansatz)
        @runtime v set_optimized!(ansatz, false)
        @runtime v set_converged!(ansatz, false)

        @runtime v angles(ansatz)
        x = zeros(typeof_parameter(ansatz), length(ansatz))
        @runtime v bind!(ansatz, x)

        @runtime v ansatz[1]
        @runtime v ansatz[1] = (generator => angle)
    end

    v && println("Runtime from `__protocols`:")
    @testset "Protocols" begin
        @runtime v typeof_score(adapt)
        @runtime v calculate_score(ansatz, adapt, generator, observable, reference)
        @runtime v calculate_scores(ansatz, adapt, pool, observable, reference)

        ansatz = deepcopy(ansatz)   # Avoid mutating the test ansatz.
        @runtime v adapt!(ansatz, Trace(),
                adapt, pool, observable, reference, callbacks)
        @runtime v optimize!(ansatz, Trace(),
                vqe, observable, reference, callbacks)
        @runtime v run!(ansatz, Trace(),
                adapt, vqe, pool, observable, reference, callbacks)
    end

    end # @testset "Runtime Validation"
end


"""
    validate_consistency(
        ansatz::AbstractAnsatz,
        adapt::AdaptProtocol,
        pool::GeneratorList,
        observable::Observable,
        reference::QuantumState,
    )

Check that every core ADAPT function is internally consistent
    (ie. different versions of the same function give consistent results).

"""
function validate_consistency(
    ansatz::AbstractAnsatz,
    adapt::AdaptProtocol,
    pool::GeneratorList,
    observable::Observable,
    reference::QuantumState,
)
    return @testset "Validate Internal Consistency" begin

    # INITIAL SETUP, to ensure we have something to test
    isempty(pool) && error("Cannot validate ADAPT with an empty pool.")
    generator = first(pool)
    angle = 1.0
    push!(ansatz, generator => angle)   # Guarantees evolution *does* something.

    @testset "Evolution by Ansatz" begin
        original = deepcopy(reference)
        evolved = evolve_state(ansatz, reference)
        @test original == reference

        evolvedinplace = evolve_state!(ansatz, original)
        @test original === evolvedinplace
        @test evolved == evolvedinplace
    end

    @testset "Evolution by Generator" begin
        original = deepcopy(reference)
        evolved = evolve_state(generator, angle, reference)
        @test original == reference

        evolvedinplace = evolve_state!(generator, angle, original)
        @test original === evolvedinplace
        @test evolved == evolvedinplace
    end

    @testset "Observable Estimation" begin
        expval = evaluate(ansatz, observable, reference)

        evolved = evolve_state(ansatz, reference)
        expval_evolved = evaluate(observable, evolved)
        @test expval == expval_evolved

        fn = make_costfunction(ansatz, observable, reference)
        x = angles(ansatz)
        @test expval == fn(x)
    end

    @testset "VQE Gradient" begin
        grad = gradient(ansatz, observable, reference)
        tobe_gradinplace = similar(grad)

        gradinplace = gradient!(tobe_gradinplace, ansatz, observable, reference)
        @test tobe_gradinplace === gradinplace
        @test grad == gradinplace

        eachpartial = [
            partial(i, ansatz, observable, reference)
                for i in eachindex(grad)
        ]
        @test grad ≈ eachpartial

        gd = make_gradfunction(ansatz, observable, reference)
        x = angles(ansatz)
        @test grad == gd(x)

        g! = make_gradfunction!(ansatz, observable, reference)
        gradinplace = g!(tobe_gradinplace, x)
        @test tobe_gradinplace === gradinplace
        @test grad == gradinplace
    end

    @testset "Ansatz Behavior" begin
        # CHECK ITERATION
        Gs = ADAPT.__get__generators(ansatz)
        θs = ADAPT.__get__parameters(ansatz)
        @test all(
            i -> ((G, θ) = ansatz[i]; (G === Gs[i]) && (θ === θs[i])),
            eachindex(ansatz),
        )

        # CHECK ANGLES
        angles = deepcopy(ADAPT.angles(ansatz))
        @test angles == ADAPT.__get__parameters(ansatz)
        newangles = .- angles
        bind!(ansatz, newangles)
        @test newangles == ADAPT.__get__parameters(ansatz)
        bind!(ansatz, angles)                   # RESET

        # CHECK `optimized` FLAG
        optimized = is_optimized(ansatz)
        @test optimized == ADAPT.__get__optimized(ansatz)[]
        notoptimized = !optimized
        set_optimized!(ansatz, notoptimized)
        @test notoptimized == ADAPT.__get__optimized(ansatz)[]
        set_optimized!(ansatz, optimized)       # RESET

        # CHECK `converged` FLAG
        converged = is_converged(ansatz)
        @test converged == ADAPT.__get__converged(ansatz)[]
        notconverged = !converged
        set_converged!(ansatz, notconverged)
        @test notconverged == ADAPT.__get__converged(ansatz)[]
        set_optimized!(ansatz, converged)       # RESET
    end

    @testset "ADAPT Scores" begin
        scores = calculate_scores(ansatz, adapt, pool, observable, reference)

        eachscore = [
            calculate_score(ansatz, adapt, pool[i], observable, reference)
                for i in eachindex(pool)
        ]
        @test scores == eachscore
    end

    end # @testset "Validate Internal Consistency"
end


"""
    validate_evolution(
        generator::Generator,
        angle::Parameter,
        reference::QuantumState;
        tolerance=1e-10,
    )

Check that generator evolution matches brute-force matrix-vector results.

The difference vector between core ADAPT and brute-force
    must have a norm within `tolerance`.

This function requires the following constructors to be defined:
- Matrix(::Generator)
- Vector(::QuantumState)

"""
function validate_evolution(
    generator::Generator,
    angle::Parameter,
    reference::QuantumState;
    tolerance=1e-10,
)
    return @testset "Validate Generator Evolution" begin

    # CORE ADAPT SOLUTION
    state = evolve_state(generator, angle, reference)

    # CONVERT DATA STRUCTURES TO ARRAYS
    ψREF = Vector(reference)
    ψTGT = Vector(state)
    G = Matrix(generator)

    # COMPUTE SOLUTION BY BRUTE-FORCE
    U = cis((-angle) .* G)
    ψ = U * ψREF

    # CHECK SOLUTIONS
    @test LinearAlgebra.norm(ψ .- ψTGT) ≤ tolerance

    end # @testset "Validate Generator Evolution"
end


"""
    validate_evaluation(
        observable::Observable,
        reference::QuantumState;
        tolerance=1e-10,
    )

Check that observable evaluation matches brute-force matrix-vector results.

The difference between core ADAPT and brute-force
    must have an absolute value within `tolerance`.

This function requires the following constructors to be defined:
- Matrix(::Observable)
- Vector(::QuantumState)

"""
function validate_evaluation(
    observable::Observable,
    reference::QuantumState;
    tolerance=1e-10,
)
    return @testset "Validate Observable Estimation" begin

    # CORE ADAPT SOLUTION
    energy = evaluate(observable, reference)

    # CONVERT DATA STRUCTURES TO ARRAYS
    ψ = Vector(reference)
    H = Matrix(observable)

    # COMPUTE SOLUTION BY BRUTE-FORCE
    E = ψ' * H * ψ

    # CHECK SOLUTIONS
    @test abs(E - energy) ≤ tolerance

    end # @testset "Validate Observable Estimation"
end

"""
    validate_gradient(
        ansatz::AbstractAnsatz,
        observable::Observable,
        reference::QuantumState;
        tolerance=1e-10,
    )

Check that the gradient function matches the finite difference.

The difference vector between core ADAPT and brute-force
    must have a norm within `tolerance`.

"""
function validate_gradient(
    ansatz::AbstractAnsatz,
    observable::Observable,
    reference::QuantumState;
    tolerance=1e-10,
)
    return @testset "Validate Gradient" begin

    # CORE ADAPT SOLUTION
    g0 = gradient(ansatz, observable, reference)

    # SETUP FOR THE FINITE DIFFERENCE
    cfd = FiniteDifferences.central_fdm(5, 1)
    f = x -> (
        xi = deepcopy(angles(ansatz));
        bind!(ansatz, x);
        E = evaluate(ansatz, observable, reference);
        bind!(ansatz, xi);
        E
    )
    x0 = angles(ansatz)

    # RUN THE FINITE DIFFERENCE
    gΔ = FiniteDifferences.grad(cfd, f, x0)[1]

    # CHECK SOLUTIONS
    @test LinearAlgebra.norm(g0 .- gΔ) ≤ tolerance

    end # @testset "Validate Gradient"
end


"""
    validate_score(
        ansatz::AbstractAnsatz,
        adapt::AdaptProtocol,
        pool::GeneratorList,
        observable::Observable,
        reference::QuantumState;
        tolerance=1e-10,
    )

Check that the score for each pool operator matches
    the partial for that pool operator when added to a candidate ansatz.

Of course this only makes sense when the score is the gradient,
    which depends on the ADAPT protocol.
But this is a common-enough choice to justify a standard method.
Other ADAPT protocols may override this method, if desired.

The difference vector between core ADAPT and brute-force
    must have a norm within `tolerance`.

"""
function validate_scores(
    ansatz::AbstractAnsatz,
    adapt::AdaptProtocol,
    pool::GeneratorList,
    observable::Observable,
    reference::QuantumState;
    tolerance=1e-10,
)
    return @testset "Validate Scores" begin

    # CORE ADAPT SOLUTION
    scores = calculate_scores(ansatz, adapt, pool, observable, reference)

    # CONSTRUCT SOLUTION FROM PARTIALS
    L = length(ansatz)
    partials = similar(scores)
    for (i, generator) in enumerate(pool)
        candidate = deepcopy(ansatz)
        push!(candidate, generator => zero(typeof_parameter(candidate)))
        partials[i] = partial(L+1, candidate, observable, reference)
    end

    # CHECK SOLUTIONS
    @test LinearAlgebra.norm(scores .- abs.(partials)) ≤ tolerance

    end # @testset "Validate Gradient"
end















#= TODO: I've decided to remove (eventually) __matrix_functions from core,
    so its tests perhaps belong elsewhere. Maybe even inside its module? =#

# function validate_runtime_matrix(
#     ansatz::AbstractAnsatz,
#     ADAPT::AdaptProtocol,
#     VQE::OptimizationProtocol,
#     pool::GeneratorList,
#     observable::Observable,
#     reference::QuantumState,
#     callbacks::CallbackList,
# )
#     # INITIAL SETUP, to ensure we have something to test
#     isempty(pool) && error("Cannot validate ADAPT with an empty pool.")
#     generator = first(pool)
#     angle = 1.0
#     push!(ansatz, generator => angle)   # Guarantees evolution *does* something.

#     #= Not exactly sure how I want to do this yet,
#             but I've decided that `__check_matrix` should NOT be a core file.
#     Rather, this is an extra module supported for quantum states
#             that support some kind of `as_vector()` function.
#     Making a unitary matrix only makes sense if it acts on a particular kind of vector,
#             but the kind of vector may depend on the type of the reference
#             (eg. a density matrix might be represented as a superket).
#     So, these methods need to be dispatched on a QuantumState!
#     (Of course they can default to assuming the QuantumState is a statevector,
#             so the methods we have at present may be fine,
#             but I'll need to rethink the interface eventually.)
#     =#
#     vector = Vector(reference)  # TODO: I suspect I'll want a new function for this.
#     N = size(vector, 1)
#     unitary = Matrix{ComplexF64}(LinearAlgebra.I, N, N)
#     # METHODS FROM `__matrix_functions`
#     @runtime v evolve_unitary(ansatz, unitary)
#     @runtime v evolve_unitary(generator, angle, unitary)
#     @runtime v evolve_unitary!(ansatz, unitary)
#     @runtime v evolve_unitary!(generator, angle, unitary)
#     @runtime v Matrix(Float64, N, ansatz)
#     @runtime v Matrix(N, ansatz)
#     @runtime v Matrix(Float64, N, generator, angle)
#     @runtime v Matrix(N, generator, angle)
# end

