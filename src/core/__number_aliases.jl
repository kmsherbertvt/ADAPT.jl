#= DISCLAIMER: This file is meant to be included in ADAPT.jl.
    Names not defined in this file or in base Julia can be found there.
=#

"""
    Parameter

Semantic alias for the coefficient of a generator in an ansatz.

"""
Parameter = Number

"""
    ParameterList

Semantic alias for a vector of parameters.

"""
ParameterList = AbstractVector{<:Parameter}

"""
    Energy

Semantic alias for the expectation value of an observable.

"""
Energy = Number

"""
    EnergyList

Semantic alias for a vector of energies.

"""
EnergyList = AbstractVector{<:Energy}

"""
    Score

Semantic alias for the importance-factor of a pool operator,
    eg. the expecation value of its commutator with an observable.

"""
Score = Number

"""
    ScoreList

Semantic alias for a vector of scores.

"""
ScoreList = AbstractVector{<:Score}
