import ..ADAPT
import PauliOperators: ScaledPauliVector

#=

TODO:
The great desire: to write a QAOA ansatz compatible with TETRIS, without rewriting TETRIS.

But...we CAN rewrite TETRIS... We don't need a new type; just override the adapt! function.
That's how multiple dispatch is supposed to work. :)



Easiest re-spec: can we have a *sentinel operator* of type G,
    but which when encountered evolves the QAOAObservable instead..?
If so, we don't need separate vectors for β and γ, though we'll want convenience methods.
This seems like the most honest way to preserve the [G => F] structure.
But what should be the sentinel? I don't actually like this.


I'd probably prefer a pure sequence of [G=>F] for just the mixers,
    but with γ_parameters and γ_indices indicating where in the sequence the QAOAObervable should be evolved.
A single update function layer!() γ0 to γ_parameters and adds length(β_parameters) to γ_indices. Adding this call is the only edit needed in the new adapt! method.

But this requires rewriting evolve_state.
Well that's probably better anyway...

=#

"""
    TetrisQAOAAnsatz{F<:Parameter,G<:Generator}(
        observable::QAOAObservable,
        γ0::F,
        generators::Vector{G},
        β_parameters::Vector{F},
        γ_parameters::Vector{F},
        optimized::Bool,
        converged::Bool,
    )

An ADAPT state suitable for ADAPT-QAOA.
The standard ADAPT generators are interspersed with the observable itself.

# Type Parameters
- `F`: the number type for the parameters (usually `Float64` is appropriate).
- `G`: the generator type.

# Parameter
- `observable`: the observable, which is interspersed with generators when evolving
- `γ0`: initial coefficient of the observable, whenever a new generator is added
- `generators`: list of current generators (i.e. mixers)
- `β_parameters`: list of current generator coefficients
- `γ_parameters`: list of current observable coefficients
- `optimized`: whether the current parameters are flagged as optimal
- `converged`: whether the current generators are flagged as converged

"""
struct TetrisQAOAAnsatz{F,G} <: ADAPT.AbstractAnsatz{F,G}
    observable::QAOAObservable
    γ0::F
    γ_values::Vector{F}
    γ_layers::Vector{Int}
    generators::Vector{G}
    parameters::Vector{F}
    optimized::Ref{Bool}
    converged::Ref{Bool}
end

function insertlayer!(ansatz::TetrisQAOAAnsatz)
    push!(ansatz.γ_values, γ0)
    push!(ansatz.γ_layers, 1+length(ansatz.parameters))
    return ansatz
end

function nlayers(ansatz::TetrisQAOAAnsatz)
    return length(ansatz.γ_values)
end

ADAPT.__get__optimized(ansatz::TetrisQAOAAnsatz) = ansatz.optimized
ADAPT.__get__converged(ansatz::TetrisQAOAAnsatz) = ansatz.converged
#= NOTE: We implement the following getters because the interface requires it,
    but we are sort of redefining in this ansatz what they mean
    (e.g. `__get__parameters` does not include all parameters, only includes β angles),
    and overriding all the methods where they are used. =#
ADAPT.__get__generators(ansatz::TetrisQAOAAnsatz) = ansatz.generators
ADAPT.__get__parameters(ansatz::TetrisQAOAAnsatz) = ansatz.parameters




function ADAPT.angles(ansatz::TetrisQAOAAnsatz{F,G}) where {F,G}
    #= TODO (but not a priority and I don't know how):
    We should be able to make this a view, to avoid allocations.
    If/when done, some `copy(angles)` should become `collect(angles)`.
    =#
    return vcat(ansatz.γ_values, ansatz.parameters)
end

function ADAPT.bind!(ansatz::TetrisQAOAAnsatz{F,G}, x::AbstractVector{F}) where {F,G}
    L = nlayers(ansatz)
    ansatz.γ_values .= @view(x[1:L])
    ansatz.parameters .= @view(x[L+1:end])
    return ansatz
end

function ADAPT.evolve_state!(ansatz::TetrisQAOAAnsatz, state::QuantumState)
    L = nlayers(ansatz)
    l = 0
    for i in eachindex(ansatz.parameters)
        while l ≤ L && ansatz.γ_layers[l] == i
            ADAPT.evolve_state!(ansatz.observable, ansatz.γ_values[l], state)
            l += 1
        end
        ADAPT.evolve_state!(ansatz.observable, ansatz.γ_values[l], state)
    end
    while l ≤ L # ansatz.γ_layers[l] > max i
        ADAPT.evolve_state!(ansatz.observable, ansatz.γ_values[l], state)
        l += 1
    end
    return state
end








##########################################################################################

function ADAPT.gradient!(
    result::AbstractVector,
    ansatz::TetrisQAOAAnsatz,
    observable::AnyPauli,
    reference::ADAPT.QuantumState,
)
    #= TODO: This function so far has just been copied from basics/pauliplugin.jl.

    I've landed on an ansatz structure where `result` is a vector to be filled
        with all γ parameter partials, followed by all β parameter partials.

    You'll probably need to copy/paste/modify the symbols `AnyPauli`
        and `__make__costate` from the same location.

    The function will need to be modified
        to account for the QAOAObservable in the reverse evolution,
        and of course the gradients for the γ parameters also need to be computed,
        perhaps using a totally different costate function. :/
    I haven't thought it through.

    I suggest starting by staring at the `evolve_state!` function above,
        to understand how forward evolution is supposed to work.
    (Maybe even better to start by testing the `evolve_state!` function
        to make sure it works as intended...)

    =#

    ψ = ADAPT.evolve_state(ansatz, reference)   # FULLY EVOLVED ANSATZ |ψ⟩
    λ = observable * ψ                          # CALCULATE |λ⟩ = H |ψ⟩

    for i in reverse(eachindex(ansatz))
        G, θ = ansatz[i]
        ADAPT.evolve_state!(G', -θ, ψ)          # UNEVOLVE BRA
        σ = __make__costate(G, θ, ψ)            # CALCULATE |σ⟩ = exp(-iθG) (-iG) |ψ⟩
        result[i] = 2 * real(dot(σ, λ))         # CALCULATE GRADIENT ⟨λ|σ⟩ + h.t.
        ADAPT.evolve_state!(G', -θ, λ)          # UNEVOLVE KET
    end

    return result
end

##########################################################################################

#= Make compatible with TETRIS protocol. =#

function ADAPT.adapt!(
    ansatz::TetrisQAOAAnsatz,
    trace::ADAPT.Trace,
    adapt_type::ADAPT.TETRIS_ADAPT.TETRISADAPT,
    pool::ADAPT.GeneratorList,
    observable::ADAPT.Observable,
    reference::ADAPT.QuantumState,
    callbacks::ADAPT.CallbackList,
)
    adapted = invoke(ADAPT.adapt!, Tuple{
        ADAPT.AbstractAnsatz,
        ADAPT.Trace,
        ADAPT.TETRIS_ADAPT.TETRISADAPT,
        ADAPT.GeneratorList,
        ADAPT.Observable,
        ADAPT.QuantumState,
        ADAPT.CallbackList,
    }, ansatz, trace, adapt_type, pool, observable, reference, callbacks)
    adapted || return adapted   # STOP IF ADAPTATION IS TERMINATED

    # ADD A NEW LAYER
    insertlayer!(ansatz)
    return adapted
end

function ADAPT.calculate_score(
    ansatz::TetrisQAOAAnsatz,
    ::ADAPT.TETRIS_ADAPT.TETRISADAPT,
    generator::AnyPauli,
    observable::AnyPauli,
    reference::ADAPT.QuantumState,
)
    state = ADAPT.evolve_state(ansatz, reference)
    ADAPT.evolve_state!(ansatz.observable, ansatz.γ0, state)
    return abs(ADAPT.Basics.MyPauliOperators.measure_commutator(
            generator, observable, state))
end
